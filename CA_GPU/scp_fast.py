"""GPU/bitset-accelerated front-end for the exact Coverage-Axis SCP.

Keeps the exact solvers from scp_exact (CP-SAT / HiGHS) and the same global
optimum, with the front-end work in bitset / GPU form:

  * coverage relation built as packed bitsets on the GPU (Warp): for each
    surface point, the set of candidate balls containing it (rowbits), and the
    transpose (colbits).  ||s_j - c_i|| < r_i, matching the reference.
  * vectorized exact presolve: unique-cover forcing (+ removing rows a forced
    candidate already covers) and identical-row/column dedup, all numpy bitset
    ops, optimum-preserving.
  * vectorized bitset greedy for the warm-start incumbent.

The heavier O(n^2)/O(m^2) dominance reductions live in scp_exact.presolve and
remain available (flag); they are off by default here, as CP-SAT's own presolve
already exploits dominance.
"""
from __future__ import annotations
import time
import numpy as np

import scp_exact as SX

try:
    import warp as wp
    _HAS_WARP = True
except Exception:
    _HAS_WARP = False


# ----------------------------------------------------------------------------
# GPU coverage -> packed bitsets
# ----------------------------------------------------------------------------
if _HAS_WARP:
    wp.set_module_options({"enable_backward": False})

    @wp.kernel
    def _rowbits_kernel(surf: wp.array2d(dtype=wp.float64),
                        cen: wp.array2d(dtype=wp.float64),
                        r2: wp.array(dtype=wp.float64),
                        n: wp.int32, Wn: wp.int32,
                        rowbits: wp.array2d(dtype=wp.uint64)):
        j = wp.tid()                       # one thread per surface point
        sx = surf[j, 0]; sy = surf[j, 1]; sz = surf[j, 2]
        for i in range(n):
            dx = sx - cen[i, 0]; dy = sy - cen[i, 1]; dz = sz - cen[i, 2]
            if dx * dx + dy * dy + dz * dz < r2[i]:
                w = i >> 6
                b = wp.uint64(1) << wp.uint64(i & 63)
                rowbits[j, w] = rowbits[j, w] | b

    @wp.kernel
    def _colbits_kernel(surf: wp.array2d(dtype=wp.float64),
                        cen: wp.array2d(dtype=wp.float64),
                        r2: wp.array(dtype=wp.float64),
                        m: wp.int32, Wm: wp.int32,
                        colbits: wp.array2d(dtype=wp.uint64)):
        i = wp.tid()                       # one thread per candidate
        cx = cen[i, 0]; cy = cen[i, 1]; cz = cen[i, 2]; ri2 = r2[i]
        for j in range(m):
            dx = surf[j, 0] - cx; dy = surf[j, 1] - cy; dz = surf[j, 2] - cz
            if dx * dx + dy * dy + dz * dz < ri2:
                w = j >> 6
                b = wp.uint64(1) << wp.uint64(j & 63)
                colbits[i, w] = colbits[i, w] | b


def warmup():
    """Trigger the one-time wp.init()/CUDA context + kernel compile so it is not
    charged to the first solve (it is amortized in any real batch use)."""
    if _HAS_WARP and wp.is_cuda_available():
        build_bitsets(np.zeros((2, 3)), np.zeros((2, 3)), np.ones(2), use_gpu=True)


def build_bitsets(surface, centers, radii, use_gpu=True):
    """Return (rowbits[m,Wn], colbits[n,Wm]) uint64 packed coverage bitsets.
    Covered iff ||s-c||^2 < r^2 (strict), matching the reference."""
    m, n = len(surface), len(centers)
    Wn, Wm = (n + 63) // 64, (m + 63) // 64
    r2 = (radii.astype(np.float64) ** 2)
    if use_gpu and _HAS_WARP and wp.is_cuda_available():
        dev = "cuda:0"
        s_d = wp.array(surface.astype(np.float64), dtype=wp.float64, device=dev)
        c_d = wp.array(centers.astype(np.float64), dtype=wp.float64, device=dev)
        r_d = wp.array(r2, dtype=wp.float64, device=dev)
        rb = wp.zeros((m, Wn), dtype=wp.uint64, device=dev)
        cb = wp.zeros((n, Wm), dtype=wp.uint64, device=dev)
        wp.launch(_rowbits_kernel, dim=m, inputs=[s_d, c_d, r_d, n, Wn, rb], device=dev)
        wp.launch(_colbits_kernel, dim=n, inputs=[s_d, c_d, r_d, m, Wm, cb], device=dev)
        wp.synchronize_device(dev)
        return rb.numpy(), cb.numpy()
    # numpy fallback (chunked)
    rowbits = np.zeros((m, Wn), dtype=np.uint64)
    colbits = np.zeros((n, Wm), dtype=np.uint64)
    for s in range(0, m, 512):
        e = min(s + 512, m)
        cov = ((surface[s:e, None, :] - centers[None, :, :]) ** 2).sum(-1) < r2[None, :]
        for k in range(e - s):
            ci = np.nonzero(cov[k])[0]
            j = s + k
            for i in ci:
                rowbits[j, i >> 6] |= np.uint64(1) << np.uint64(int(i) & 63)
                colbits[i, j >> 6] |= np.uint64(1) << np.uint64(j & 63)
    return rowbits, colbits


# ----------------------------------------------------------------------------
# vectorized helpers
# ----------------------------------------------------------------------------
def _popcount_rows(bits: np.ndarray) -> np.ndarray:
    """Per-row popcount of a [k, W] uint64 bitset array."""
    v = bits.copy()
    v = v - ((v >> np.uint64(1)) & np.uint64(0x5555555555555555))
    v = (v & np.uint64(0x3333333333333333)) + ((v >> np.uint64(2)) & np.uint64(0x3333333333333333))
    v = (v + (v >> np.uint64(4))) & np.uint64(0x0f0f0f0f0f0f0f0f)
    pc = ((v * np.uint64(0x0101010101010101)) >> np.uint64(56)).astype(np.int64)
    return pc.sum(axis=1)


def _bits_to_ids(word_row: np.ndarray) -> np.ndarray:
    """Indices of set bits in a [W] uint64 row."""
    ids = []
    for w, val in enumerate(word_row):
        x = int(val)
        base = w * 64
        while x:
            b = (x & -x).bit_length() - 1
            ids.append(base + b)
            x &= x - 1
    return np.array(ids, dtype=np.int64)


# ----------------------------------------------------------------------------
# vectorized exact presolve (optimum-preserving): unique-cover + dedup
# ----------------------------------------------------------------------------
def presolve_fast(rowbits: np.ndarray, n_cols: int):
    """Unique-cover forcing + remove satisfied rows + dedup identical rows.
    Returns (forced ids, kept_rows mask, rowbits) — all optimum-preserving."""
    m, Wn = rowbits.shape
    rb = rowbits.copy()
    # Infeasibility guard: a surface row with no covering candidate cannot be
    # satisfied -> the SCP is infeasible (matches scp_exact.presolve). Cannot
    # happen on real Coverage-Axis instances (dilation covers every sample),
    # but guard so an infeasible input never returns a bogus "optimal" cover.
    empty = (_popcount_rows(rb) == 0)
    if empty.any():
        raise ValueError(f"surface row {int(np.nonzero(empty)[0][0])} is uncovered "
                         f"-> SCP infeasible")
    alive = np.ones(m, dtype=bool)
    forced_words = np.zeros(Wn, dtype=np.uint64)
    forced_any = False
    while True:
        pc = _popcount_rows(rb)
        pc[~alive] = 0
        uniq = np.nonzero((pc == 1) & alive)[0]
        if uniq.size == 0:
            break
        # OR all unique-cover rows into the forced set
        new_forced = np.zeros(Wn, dtype=np.uint64)
        for j in uniq:
            new_forced |= rb[j]
        before = forced_words.copy()
        forced_words |= new_forced
        forced_any = True
        # remove rows that share any candidate with the forced set
        inter = (rb & forced_words[None, :])
        covered = (inter != 0).any(axis=1)
        alive &= ~covered
        if np.array_equal(before, forced_words) and not covered.any():
            break
    forced = _bits_to_ids(forced_words)
    # dedup identical alive rows (same bitset) -> keep one
    if alive.any():
        live_idx = np.nonzero(alive)[0]
        seen = {}
        for j in live_idx:
            key = rb[j].tobytes()
            if key in seen:
                alive[j] = False
            else:
                seen[key] = j
    return forced, alive, rb


# ----------------------------------------------------------------------------
# vectorized bitset greedy (warm start)
# ----------------------------------------------------------------------------
def greedy_bitset(colbits: np.ndarray, m: int, active_cols: np.ndarray | None = None):
    """Max-coverage greedy using column bitsets over m rows. Returns chosen
    ORIGINAL candidate ids (a feasible cover of all m rows)."""
    n, Wm = colbits.shape
    cols = active_cols if active_cols is not None else np.arange(n)
    covered = np.zeros(Wm, dtype=np.uint64)
    target_pc = m
    chosen = []
    cur = colbits[cols].copy()
    while True:
        # gain_i = popcount(cur_i & ~covered)
        notc = ~covered
        gain = _popcount_rows(cur & notc[None, :])
        bi = int(np.argmax(gain))
        if gain[bi] == 0:
            break
        chosen.append(int(cols[bi]))
        covered |= colbits[cols[bi]]
        if int(_popcount_rows(covered[None, :])[0]) >= target_pc:
            break
    return chosen


# ----------------------------------------------------------------------------
# pipeline
# ----------------------------------------------------------------------------
def coverage_axis_fast(surface, centers, radii, *, backend="cpsat", use_gpu=True,
                       warm=True, workers=8, time_limit=None, verbose=True):
    """GPU-accelerated exact Coverage-Axis SCP. Same global optimum as the
    reference, reached faster. Returns dict (selected original ids, objective,
    proven_optimal, timing)."""
    t0 = time.perf_counter()
    rowbits, colbits = build_bitsets(surface, centers, radii, use_gpu=use_gpu)
    m = len(surface); n = len(centers)
    t_build = time.perf_counter() - t0

    t1 = time.perf_counter()
    forced, alive, rb = presolve_fast(rowbits, n)
    fset = set(forced.tolist())
    # build reduced rows (alive rows minus forced candidates) in local col space
    live = np.nonzero(alive)[0]
    appearing = set()
    raw_rows = []
    for j in live:
        ids = [i for i in _bits_to_ids(rb[j]) if i not in fset]
        if ids:
            raw_rows.append(ids); appearing.update(ids)
    active = np.array(sorted(appearing), dtype=np.int64)
    loc = {int(c): k for k, c in enumerate(active)}
    rows_local = [[loc[i] for i in r] for r in raw_rows]
    n_local = int(active.size)
    t_presolve = time.perf_counter() - t1

    t2 = time.perf_counter()
    warmL = None
    if warm and n_local:
        # greedy over the reduced program (local), warm-start the solver
        sub_colbits = np.zeros((n_local, (len(rows_local) + 63) // 64), dtype=np.uint64)
        for jj, r in enumerate(rows_local):
            for i in r:
                sub_colbits[i, jj >> 6] |= np.uint64(1) << np.uint64(jj & 63)
        warmL = greedy_bitset(sub_colbits, len(rows_local))
    t_greedy = time.perf_counter() - t2

    t3 = time.perf_counter()
    if n_local == 0:
        sel_local, proven = [], True
    elif backend == "cpsat":
        sel_local, _, proven = SX.solve_cpsat(rows_local, n_local, warm=warmL,
                                              workers=workers, time_limit=time_limit)
    elif backend == "highs":
        sel_local, _, proven = SX.solve_highs(rows_local, n_local, time_limit=time_limit)
    else:
        raise ValueError(backend)
    t_solve = time.perf_counter() - t3

    selected = np.unique(np.concatenate([
        forced, active[np.array(sel_local, np.int64)] if sel_local else np.zeros(0, np.int64)]).astype(np.int64))
    res = {"selected": selected, "objective": int(selected.size), "proven_optimal": bool(proven),
           "n_forced": int(forced.size), "n_reduced_rows": len(rows_local), "n_reduced_cols": n_local,
           "t_build_s": t_build, "t_presolve_s": t_presolve, "t_greedy_s": t_greedy,
           "t_solve_s": t_solve, "t_total_s": time.perf_counter() - t0, "backend": backend}
    if verbose:
        print(f"[fast] build(gpu={use_gpu and _HAS_WARP}) {t_build:.2f}s | presolve {t_presolve:.2f}s "
              f"(forced {forced.size}, rows {m}->{len(rows_local)}, cols {n}->{n_local}) | "
              f"greedy {t_greedy:.2f}s | solve({backend}) {t_solve:.2f}s | "
              f"obj={res['objective']} proven={proven} | total {res['t_total_s']:.2f}s")
    return res

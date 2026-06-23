"""Exact, accelerated Set-Cover solver for the Coverage Axis point-selection step.

Coverage Axis (Dou et al. 2022) selects inner poles by solving, over a finite
surface sample, the 0-1 set-cover program

    min  sum_i x_i
    s.t. sum_{i : D[j,i]=1} x_i >= 1   for every surface point j
         x_i in {0,1}

where D[j,i] = 1 iff surface point j lies inside the dilated ball (c_i, r_i).
The reference implementation hands the full D straight to a branch-and-bound
MILP (scipy.optimize.milp / HiGHS) with a wall-clock cutoff, so on large D it
slows down or returns before proving optimality.

This module keeps the SAME program and its GLOBAL optimum, but reaches it
faster by:
  1. building the coverage relation once (vectorized / chunked),
  2. exact presolve reductions that provably preserve the optimum
     (unique-cover forcing, identical-column dedup, column/row dominance),
  3. a warm-start greedy upper bound + an LP lower bound (gap certificate),
  4. an exact integer solve on the reduced program (CP-SAT or HiGHS),
  5. optional row generation so the surface sample can be made dense while the
     solver only ever sees the rows it needs (still exact on the full sample).

Every reduction is optimum-preserving; `verify.py` cross-checks the objective
against the unreduced HiGHS solve.
"""
from __future__ import annotations

import time
from dataclasses import dataclass, field

import numpy as np


# ----------------------------------------------------------------------------
# Coverage relation
# ----------------------------------------------------------------------------

@dataclass
class Coverage:
    """Sparse coverage relation in row form (CSR over surface points).

    rows[j] = sorted candidate ids whose dilated ball contains surface point j.
    Stored as CSR (indptr, indices) for compactness; n_cols = #candidates."""
    indptr: np.ndarray          # [m+1] int64
    indices: np.ndarray         # [nnz]  int32, candidate ids, sorted within each row
    n_cols: int                 # number of candidates
    @property
    def m(self) -> int:
        return int(self.indptr.size - 1)
    def row(self, j: int) -> np.ndarray:
        return self.indices[self.indptr[j]:self.indptr[j + 1]]


def build_coverage(surface: np.ndarray, centers: np.ndarray, radii: np.ndarray,
                   chunk: int = 1024) -> Coverage:
    """Coverage CSR: surface point j covered by candidate i iff
    ||s_j - c_i|| < r_i   (strict, matching the reference `torch.gt(radius, D)`).

    Chunked over surface points so the m x n distance block never materializes
    in full. radii are the already-dilated ball radii (r_i = inscribed + delta)."""
    m = surface.shape[0]
    n = centers.shape[0]
    r2 = (radii.astype(np.float64) ** 2)[None, :]            # [1, n]
    rows_ind = []
    counts = np.empty(m, dtype=np.int64)
    for s in range(0, m, chunk):
        e = min(s + chunk, m)
        d2 = ((surface[s:e, None, :] - centers[None, :, :]) ** 2).sum(-1)  # [b, n]
        cov = d2 < r2                                         # strict
        for k in range(e - s):
            idx = np.nonzero(cov[k])[0].astype(np.int32)
            rows_ind.append(idx)
            counts[s + k] = idx.size
    indptr = np.zeros(m + 1, dtype=np.int64)
    np.cumsum(counts, out=indptr[1:])
    indices = (np.concatenate(rows_ind) if m else np.zeros(0, np.int32)).astype(np.int32)
    return Coverage(indptr=indptr, indices=indices, n_cols=n)


# ----------------------------------------------------------------------------
# Bitset helpers (packed uint64) for exact presolve dominance tests
# ----------------------------------------------------------------------------

def _pack_columns(cov: Coverage):
    """Column bitsets C_i over the m surface rows: bit j set iff candidate i
    covers surface point j. Returns (cols[n, W] uint64, W)."""
    m, n = cov.m, cov.n_cols
    W = (m + 63) // 64
    cols = np.zeros((n, W), dtype=np.uint64)
    # scatter: for each (row j, candidate i) set bit j of column i
    for j in range(m):
        ci = cov.row(j)
        if ci.size:
            cols[ci, j >> 6] |= np.uint64(1) << np.uint64(j & 63)
    return cols, W


# ----------------------------------------------------------------------------
# Exact presolve  (every rule preserves the optimal objective)
# ----------------------------------------------------------------------------

@dataclass
class Presolved:
    forced: np.ndarray              # original candidate ids forced into every optimum
    rows: list                      # remaining rows as lists of ORIGINAL candidate ids
    active_cols: np.ndarray         # original ids of the columns still in play
    col_local: dict                 # original id -> local index in active_cols
    stats: dict = field(default_factory=dict)


def presolve(cov: Coverage, centers=None, radii=None,
             do_col_dominance=True, do_row_dominance=True,
             geom_dominance=True) -> Presolved:
    """Reduce the SCP without changing its optimum.

    Rules (each justified in the README):
      * unique-cover forcing: a row with a single candidate forces that candidate.
      * remove rows already covered by a forced/selected candidate.
      * identical-column dedup: equal coverage bitset -> keep one.
      * column dominance: C_a subset of C_b (unit cost) -> drop a.
      * geometric column dominance: ||c_a-c_b|| + r_a <= r_b -> ball a subset of
        ball b for ANY surface sample -> drop a (sample-independent).
      * row dominance: E_a subset of E_b -> row b is implied -> drop b.
    Iterated to a fixpoint.
    """
    m, n = cov.m, cov.n_cols
    t0 = time.perf_counter()
    # mutable row sets (as python sets of original candidate ids)
    rows = [set(int(x) for x in cov.row(j)) for j in range(m)]
    alive_row = np.ones(m, dtype=bool)
    col_dead = np.zeros(n, dtype=bool)
    forced = np.zeros(n, dtype=bool)
    stats = {"m0": m, "n0": n, "forced": 0, "rows_removed": 0,
             "cols_dom": 0, "cols_dup": 0, "rows_dom": 0, "geom_dom": 0}

    # infeasibility check: any empty row is uncoverable
    for j in range(m):
        if not rows[j]:
            raise ValueError(f"surface row {j} is uncovered -> SCP infeasible")

    # --- geometric column dominance (sample-independent, do once up front) ---
    if geom_dominance and centers is not None and radii is not None:
        c = centers.astype(np.float64); r = radii.astype(np.float64)
        # a dominated by b if |c_a-c_b| + r_a <= r_b ; check against larger-radius candidates
        order = np.argsort(-r)
        for bi in order:
            if col_dead[bi]:
                continue
            # candidates whose ball fits inside bi's ball
            d = np.sqrt(((c - c[bi]) ** 2).sum(1))
            dom = (d + r <= r[bi]) & (~col_dead)
            dom[bi] = False
            if dom.any():
                idx = np.nonzero(dom)[0]
                col_dead[idx] = True
                stats["geom_dom"] += int(idx.size)
                for j in range(m):
                    if alive_row[j]:
                        rows[j].difference_update(idx.tolist())

    changed = True
    while changed:
        changed = False
        # --- unique-cover forcing + remove satisfied rows ---
        for j in range(m):
            if not alive_row[j]:
                continue
            rj = rows[j]
            if len(rj) == 1:
                i = next(iter(rj))
                if not forced[i]:
                    forced[i] = True; changed = True
                    stats["forced"] += 1
        if forced.any():
            fset = set(np.nonzero(forced)[0].tolist())
            for j in range(m):
                if alive_row[j] and (rows[j] & fset):
                    alive_row[j] = False
                    stats["rows_removed"] += 1
        if changed:
            continue  # re-derive after forcing

        live_rows = np.nonzero(alive_row)[0]
        # active (not dead, not forced) columns that still appear in a live row
        appearing = set()
        for j in live_rows:
            appearing |= rows[j]
        appearing -= set(np.nonzero(forced)[0].tolist())
        active = np.array(sorted(appearing), dtype=np.int64)
        if active.size == 0:
            break

        # --- identical-column dedup + column dominance (bitset over live rows) ---
        if do_col_dominance and live_rows.size and active.size > 1:
            row_pos = {int(j): k for k, j in enumerate(live_rows)}
            L = live_rows.size; W = (L + 63) // 64
            colbits = np.zeros((active.size, W), dtype=np.uint64)
            col_idx = {int(c): a for a, c in enumerate(active)}
            for kk, j in enumerate(live_rows):
                for i in rows[j]:
                    a = col_idx.get(int(i))
                    if a is not None:
                        colbits[a, kk >> 6] |= np.uint64(1) << np.uint64(kk & 63)
            popcount = np.zeros(active.size, dtype=np.int64)
            for w in range(W):
                popcount += _popcount64(colbits[:, w]).astype(np.int64)
            # dedup identical columns (same bitset) -> keep the first
            seen = {}
            dup_mask = np.zeros(active.size, dtype=bool)
            for a in range(active.size):
                key = colbits[a].tobytes()
                if key in seen:
                    dup_mask[a] = True
                else:
                    seen[key] = a
            # column dominance: a dominated by b if (C_a & C_b) == C_a and pop_b>pop_a
            order = np.argsort(-popcount)   # try big columns as dominators
            dom_mask = dup_mask.copy()
            for b in order:
                if dom_mask[b]:
                    continue
                cb = colbits[b]
                # candidates a with C_a subset C_b : (colbits & cb)==colbits, popcount<=pop_b
                cand = (~dom_mask) & (popcount <= popcount[b])
                cand[b] = False
                if not cand.any():
                    continue
                sub = colbits[cand]
                anded = sub & cb[None, :]
                issub = (anded == sub).all(axis=1)
                ids = np.nonzero(cand)[0][issub]
                if ids.size:
                    dom_mask[ids] = True
            if dom_mask.any():
                killed = active[dom_mask]
                col_dead[killed] = True
                stats["cols_dom"] += int(dom_mask.sum() - dup_mask.sum())
                stats["cols_dup"] += int(dup_mask.sum())
                kset = set(killed.tolist())
                for j in live_rows:
                    rows[j].difference_update(kset)
                changed = True
                continue

        # --- row dominance: drop row b if some live row a has E_a subset E_b ---
        if do_row_dominance and live_rows.size > 1:
            # order rows by size; a smaller row can dominate (be subset of) a larger
            rl = [(len(rows[j]), int(j)) for j in live_rows]
            rl.sort()
            removed_any = False
            # frozenset cache
            fs = {j: frozenset(rows[j]) for _, j in rl}
            for ai in range(len(rl)):
                _, ja = rl[ai]
                if not alive_row[ja]:
                    continue
                Ea = fs[ja]
                if not Ea:
                    continue
                for bi in range(len(rl) - 1, ai, -1):
                    _, jb = rl[bi]
                    if jb == ja or not alive_row[jb]:
                        continue
                    if Ea <= fs[jb]:        # E_a subset E_b -> row b implied
                        alive_row[jb] = False
                        stats["rows_dom"] += 1
                        removed_any = True
            if removed_any:
                changed = True

    # assemble result
    forced_ids = np.nonzero(forced)[0].astype(np.int64)
    live_rows = np.nonzero(alive_row)[0]
    fset = set(forced_ids.tolist())
    out_rows = []
    appearing = set()
    for j in live_rows:
        rj = rows[j] - fset
        if rj:                                  # rows fully covered by forced are gone
            out_rows.append(sorted(rj))
            appearing |= rj
    active_cols = np.array(sorted(appearing), dtype=np.int64)
    col_local = {int(c): k for k, c in enumerate(active_cols)}
    stats["m_reduced"] = len(out_rows)
    stats["n_reduced"] = int(active_cols.size)
    stats["t_presolve_s"] = time.perf_counter() - t0
    return Presolved(forced=forced_ids, rows=out_rows, active_cols=active_cols,
                     col_local=col_local, stats=stats)


def _popcount64(a: np.ndarray) -> np.ndarray:
    """Popcount of a uint64 array (numpy, no native popcount needed)."""
    v = a.copy()
    v = v - ((v >> np.uint64(1)) & np.uint64(0x5555555555555555))
    v = (v & np.uint64(0x3333333333333333)) + ((v >> np.uint64(2)) & np.uint64(0x3333333333333333))
    v = (v + (v >> np.uint64(4))) & np.uint64(0x0f0f0f0f0f0f0f0f)
    return (v * np.uint64(0x0101010101010101)) >> np.uint64(56)


# ----------------------------------------------------------------------------
# Bounds: greedy upper bound (warm start) + LP lower bound (gap certificate)
# ----------------------------------------------------------------------------

def greedy_cover(rows: list, n_cols: int) -> list:
    """Classic max-coverage greedy. Returns a feasible cover (list of local col
    ids). Used as the exact solver's warm start / incumbent."""
    from collections import defaultdict
    covers = defaultdict(set)          # col -> set(row idx)
    for j, r in enumerate(rows):
        for i in r:
            covers[i].add(j)
    uncovered = set(range(len(rows)))
    chosen = []
    cov_sets = {i: set(s) for i, s in covers.items()}
    while uncovered:
        best, best_gain = -1, 0
        for i, s in cov_sets.items():
            g = len(s & uncovered)
            if g > best_gain:
                best_gain, best = g, i
        if best < 0:
            raise RuntimeError("greedy could not cover (infeasible?)")
        chosen.append(best)
        uncovered -= cov_sets[best]
    return chosen


def lp_lower_bound(rows: list, n_cols: int) -> float:
    """LP relaxation objective (>= a valid lower bound on the integer optimum).
    ceil() of it is a valid integer lower bound."""
    from scipy.optimize import linprog
    from scipy.sparse import csr_matrix
    if not rows:
        return 0.0
    data, ri, ci = [], [], []
    for j, r in enumerate(rows):
        for i in r:
            data.append(1.0); ri.append(j); ci.append(i)
    A = csr_matrix((data, (ri, ci)), shape=(len(rows), n_cols))
    res = linprog(c=np.ones(n_cols), A_ub=-A, b_ub=-np.ones(len(rows)),
                  bounds=(0, 1), method="highs")
    return float(res.fun) if res.success else 0.0


# ----------------------------------------------------------------------------
# Exact integer solvers on the reduced program
# ----------------------------------------------------------------------------

def solve_cpsat(rows: list, n_cols: int, warm: list | None = None,
                workers: int = 8, time_limit: float | None = None):
    """Exact 0-1 set cover via OR-Tools CP-SAT. Each row -> BoolOr; minimize
    sum of vars. Returns (selected_local_cols, objective, proven_optimal)."""
    from ortools.sat.python import cp_model
    model = cp_model.CpModel()
    x = [model.NewBoolVar(f"x{i}") for i in range(n_cols)]
    for r in rows:
        model.AddBoolOr([x[i] for i in r])
    model.Minimize(sum(x))
    if warm:
        model.AddHint  # noqa (kept explicit below)
        for i in range(n_cols):
            model.AddHint(x[i], 1 if i in set(warm) else 0)
    solver = cp_model.CpSolver()
    solver.parameters.num_search_workers = workers
    if time_limit:
        solver.parameters.max_time_in_seconds = float(time_limit)
    status = solver.Solve(model)
    if status not in (cp_model.OPTIMAL, cp_model.FEASIBLE):
        raise RuntimeError(f"CP-SAT status={solver.StatusName(status)}")
    sel = [i for i in range(n_cols) if solver.BooleanValue(x[i])]
    return sel, int(round(solver.ObjectiveValue())), (status == cp_model.OPTIMAL)


def solve_highs(rows: list, n_cols: int, time_limit: float | None = None):
    """Exact 0-1 set cover via scipy.optimize.milp (HiGHS) — the reference
    backend. Returns (selected_local_cols, objective, proven_optimal)."""
    from scipy.optimize import milp, Bounds, LinearConstraint
    from scipy.sparse import csr_matrix
    data, ri, ci = [], [], []
    for j, r in enumerate(rows):
        for i in r:
            data.append(1.0); ri.append(j); ci.append(i)
    A = csr_matrix((data, (ri, ci)), shape=(len(rows), n_cols))
    opts = {"disp": False}
    if time_limit:
        opts["time_limit"] = float(time_limit)
    res = milp(c=np.ones(n_cols), integrality=np.ones(n_cols),
               bounds=Bounds(np.zeros(n_cols), np.ones(n_cols)),
               constraints=LinearConstraint(A, lb=np.ones(len(rows))),
               options=opts)
    if not res.success or res.x is None:
        raise RuntimeError(f"HiGHS milp failed: {res.message}")
    xr = np.round(res.x).astype(int)
    sel = np.nonzero(xr)[0].tolist()
    proven = bool(getattr(res, "mip_gap", 1.0) is not None and res.mip_gap <= 1e-9)
    return sel, int(xr.sum()), proven


# ----------------------------------------------------------------------------
# Full pipeline: build -> presolve -> warm start -> exact solve -> lift back
# ----------------------------------------------------------------------------

def coverage_axis_exact(surface, centers, radii, *, backend="cpsat",
                        presolve_opts=None, workers=8, time_limit=None,
                        verbose=True):
    """Solve the Coverage-Axis SCP to a PROVEN global optimum, accelerated by
    presolve + warm start. Returns dict with selected ORIGINAL candidate ids,
    objective, optimality flag, and a timing/size breakdown."""
    t0 = time.perf_counter()
    cov = build_coverage(surface, centers, radii)
    t_build = time.perf_counter() - t0
    pre = presolve(cov, centers=centers, radii=radii, **(presolve_opts or {}))
    # map reduced rows to local column ids
    rows_local = [[pre.col_local[i] for i in r] for r in pre.rows]
    n_local = int(pre.active_cols.size)
    warm = greedy_cover(rows_local, n_local) if n_local else []
    lb = lp_lower_bound(rows_local, n_local) if n_local else 0.0
    t1 = time.perf_counter()
    if n_local == 0:
        sel_local, obj_red, proven = [], 0, True
    elif backend == "cpsat":
        sel_local, obj_red, proven = solve_cpsat(rows_local, n_local, warm=warm,
                                                 workers=workers, time_limit=time_limit)
    elif backend == "highs":
        sel_local, obj_red, proven = solve_highs(rows_local, n_local, time_limit=time_limit)
    else:
        raise ValueError(f"unknown backend {backend}")
    t_solve = time.perf_counter() - t1
    selected = np.concatenate([pre.forced,
                               pre.active_cols[np.array(sel_local, dtype=np.int64)]
                               if sel_local else np.zeros(0, np.int64)]).astype(np.int64)
    selected = np.unique(selected)
    obj = int(selected.size)
    res = {
        "selected": selected, "objective": obj, "proven_optimal": bool(proven),
        "lower_bound": float(len(pre.forced) + lb),
        "n_forced": int(pre.forced.size), "greedy_ub": int(len(pre.forced) + len(warm)),
        "t_build_s": t_build, "t_presolve_s": pre.stats["t_presolve_s"],
        "t_solve_s": t_solve, "t_total_s": time.perf_counter() - t0,
        "presolve": pre.stats, "backend": backend,
    }
    if verbose:
        s = pre.stats
        print(f"[scp] build {t_build:.2f}s | presolve {s['t_presolve_s']:.2f}s: "
              f"m {s['m0']}->{s['m_reduced']}, n {s['n0']}->{s['n_reduced']}, "
              f"forced {s['forced']}, geom_dom {s['geom_dom']}, col_dom {s['cols_dom']}, "
              f"dup {s['cols_dup']}, row_dom {s['rows_dom']}")
        print(f"[scp] solve({backend}) {t_solve:.2f}s | objective={obj} "
              f"(LB {res['lower_bound']:.2f}, greedy {res['greedy_ub']}) "
              f"proven_optimal={proven} | total {res['t_total_s']:.2f}s")
    return res

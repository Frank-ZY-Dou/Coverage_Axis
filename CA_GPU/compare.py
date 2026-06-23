"""Head-to-head comparison for the Coverage-Axis SCP: solution + runtime.

Separates the two costs (build vs solve) for a like-for-like comparison:

  COVERAGE BUILD   numpy-chunked        vs   GPU/Warp bitset (warmed)
  EXACT SOLVE      HiGHS milp (ref)     vs   CP-SAT          vs  presolve+CP-SAT

Every solve is exact; the objective is identical across all of them
(acceleration must not change the answer). An independent feasibility re-check
confirms each returned selection covers all surface points.
"""
from __future__ import annotations
import argparse, time
import numpy as np
import scp_exact as SX
import scp_fast as SF


def feasible(surface, centers, radii, sel):
    if sel.size == 0:
        return len(surface) == 0
    ok = np.zeros(len(surface), dtype=bool)
    cen = centers[sel]; r2 = radii[sel] ** 2
    for s in range(0, len(surface), 512):
        e = min(s + 512, len(surface))
        d2 = ((surface[s:e, None, :] - cen[None, :, :]) ** 2).sum(-1)
        ok[s:e] = (d2 < r2[None, :]).any(1)
    return bool(ok.all())


def rows_from_bitsets(rowbits, n):
    return [SF._bits_to_ids(rowbits[j]).tolist() for j in range(rowbits.shape[0])]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--inst", required=True)
    ap.add_argument("--time-limit", type=float, default=None)
    ap.add_argument("--workers", type=int, default=8)
    args = ap.parse_args()
    inst = np.load(args.inst)
    surf, cen, rad = inst["surface"], inst["centers"], inst["radii"]
    print(f"=== instance: surface={len(surf)} candidates={len(cen)} ===\n")

    # ---------- COVERAGE BUILD: numpy vs GPU ----------
    t = time.perf_counter(); cov = SX.build_coverage(surf, cen, rad); t_np = time.perf_counter() - t
    SF.warmup()                                            # amortize one-time wp.init
    t = time.perf_counter(); rb, cb = SF.build_bitsets(surf, cen, rad, use_gpu=True); t_gpu = time.perf_counter() - t
    rows = [list(int(i) for i in cov.row(j)) for j in range(cov.m)]
    rows_gpu = rows_from_bitsets(rb, cov.n_cols)
    build_match = all(set(a) == set(b) for a, b in zip(rows, rows_gpu))
    print("--- COVERAGE BUILD ---")
    print(f"  numpy-chunked : {t_np:6.2f}s")
    print(f"  GPU/Warp      : {t_gpu:6.2f}s   ({t_np/max(t_gpu,1e-9):.1f}x)   matches numpy: {build_match}\n")

    sols = {}
    # ---------- EXACT SOLVE on the same coverage ----------
    t = time.perf_counter(); selR, objR, prR = SX.solve_highs(rows, cov.n_cols, time_limit=args.time_limit); tR = time.perf_counter() - t
    sols["HiGHS milp (reference)"] = (objR, prR, tR, np.array(selR, np.int64))

    t = time.perf_counter(); selC, objC, prC = SX.solve_cpsat(rows, cov.n_cols, workers=args.workers, time_limit=args.time_limit); tC = time.perf_counter() - t
    sols["CP-SAT"] = (objC, prC, tC, np.array(selC, np.int64))

    # presolve(fast) + CP-SAT
    rF = SF.coverage_axis_fast(surf, cen, rad, backend="cpsat", use_gpu=True,
                               workers=args.workers, time_limit=args.time_limit, verbose=False)
    sols["presolve(fast)+CP-SAT"] = (rF["objective"], rF["proven_optimal"],
                                     rF["t_presolve_s"] + rF["t_greedy_s"] + rF["t_solve_s"], rF["selected"])

    print("--- EXACT SOLVE (objective = #poles; identical => acceleration is exact) ---")
    print(f"  {'method':28s} {'obj':>5} {'proven':>7} {'feasible':>9} {'solve(s)':>9} {'vs HiGHS':>9}")
    objs = []
    for name, (obj, pr, tt, sel) in sols.items():
        feas = feasible(surf, cen, rad, sel); objs.append(obj)
        print(f"  {name:28s} {obj:>5} {str(pr):>7} {str(feas):>9} {tt:>9.2f} {tR/max(tt,1e-9):>8.1f}x")
    same = len(set(objs)) == 1
    print(f"  -> all objectives identical: {'YES (exact, proven-optimal)' if same else 'NO <-- MISMATCH'}")

    print("\n--- END-TO-END (GPU build + best exact solve) vs reference (numpy build + HiGHS) ---")
    ref_e2e = t_np + tR
    fast_e2e = t_gpu + sols["CP-SAT"][2]
    print(f"  reference : {ref_e2e:6.2f}s  (numpy {t_np:.2f} + HiGHS {tR:.2f})")
    print(f"  accel     : {fast_e2e:6.2f}s  (GPU {t_gpu:.2f} + CP-SAT {sols['CP-SAT'][2]:.2f})   "
          f"=> {ref_e2e/max(fast_e2e,1e-9):.1f}x, same optimum")
    if not same:
        raise SystemExit(1)


if __name__ == "__main__":
    main()

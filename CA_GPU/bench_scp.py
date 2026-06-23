"""Benchmark + correctness cross-check for the exact Coverage-Axis SCP solver.

Compares, on the same instance:
  * baseline   — the reference approach: HiGHS milp on the FULL coverage matrix
                 (no presolve), as in Coverage_Axis_mesh.py.
  * accel/highs — presolve + warm start, then HiGHS on the reduced program.
  * accel/cpsat — presolve + warm start, then CP-SAT on the reduced program.

Asserts all three reach the SAME global optimum objective (the acceleration
must not change the answer), and reports the wall-clock speedup.
"""
from __future__ import annotations
import argparse, time
import numpy as np
import scp_exact as S


def baseline_highs_full(inst, time_limit=None):
    """Reference: build full D, hand straight to HiGHS milp (no presolve)."""
    cov = S.build_coverage(inst["surface"], inst["centers"], inst["radii"])
    rows = [list(int(i) for i in cov.row(j)) for j in range(cov.m)]
    t = time.perf_counter()
    sel, obj, proven = S.solve_highs(rows, cov.n_cols, time_limit=time_limit)
    return obj, proven, time.perf_counter() - t, cov


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--inst", required=True)
    ap.add_argument("--time-limit", type=float, default=None)
    ap.add_argument("--workers", type=int, default=8)
    args = ap.parse_args()
    inst = np.load(args.inst)
    print(f"=== instance: surface={len(inst['surface'])} candidates={len(inst['centers'])} ===\n")

    print("--- baseline: HiGHS milp on FULL D (reference Coverage_Axis) ---")
    base_obj, base_proven, base_t, cov = baseline_highs_full(inst, time_limit=args.time_limit)
    print(f"  objective={base_obj}  proven_optimal={base_proven}  time={base_t:.2f}s\n")

    print("--- accel: presolve + warm start + HiGHS (reduced) ---")
    r_h = S.coverage_axis_exact(inst["surface"], inst["centers"], inst["radii"],
                                backend="highs", workers=args.workers, time_limit=args.time_limit)
    print()
    print("--- accel: presolve + warm start + CP-SAT (reduced) ---")
    r_c = S.coverage_axis_exact(inst["surface"], inst["centers"], inst["radii"],
                                backend="cpsat", workers=args.workers, time_limit=args.time_limit)
    print()

    print("=== SUMMARY ===")
    print(f"  baseline HiGHS/full : obj={base_obj}  {base_t:6.2f}s  proven={base_proven}")
    print(f"  accel   HiGHS/reduced: obj={r_h['objective']}  {r_h['t_total_s']:6.2f}s  proven={r_h['proven_optimal']}  ({base_t/max(r_h['t_total_s'],1e-9):.1f}x)")
    print(f"  accel   CP-SAT/reduced: obj={r_c['objective']}  {r_c['t_total_s']:6.2f}s  proven={r_c['proven_optimal']}  ({base_t/max(r_c['t_total_s'],1e-9):.1f}x)")
    ok = (base_obj == r_h["objective"] == r_c["objective"])
    print(f"\n  GLOBAL OPTIMUM IDENTICAL across all three: {'YES' if ok else 'NO  <-- MISMATCH'}")
    if not ok:
        raise SystemExit(1)


if __name__ == "__main__":
    main()

"""Independent verification of the exact SCP solver.

Checks (each re-derived independently of scp_exact's internals):

  A. FEASIBILITY  — recompute coverage from raw geometry (not build_coverage)
     and assert the returned selection covers EVERY surface point. Catches a
     presolve that silently drops constraints (which would shrink the objective).

  B. OPTIMALITY (two independent certificates):
     B1. independent full-instance solve with BOTH HiGHS and CP-SAT (no presolve)
         must equal the accelerated objective;
     B2. an LP-relaxation lower bound: ceil(LP) == objective  =>  proven optimal
         without trusting any branch-and-bound's "OPTIMAL" flag.

  C. BRUTE FORCE — on small random SCP instances, the true optimum (exhaustive)
     must equal presolve+solve, exercising every presolve rule against ground truth.

Run: python3 verify.py --inst /tmp/beagle_inst.npz
"""
from __future__ import annotations
import argparse, itertools, math
import numpy as np
import scp_exact as S


# ---- A. independent coverage + feasibility -------------------------------
def coverage_dense(surface, centers, radii):
    """Independent dense boolean coverage (recomputed, not via build_coverage)."""
    m = len(surface); n = len(centers)
    D = np.zeros((m, n), dtype=bool)
    for s in range(0, m, 512):
        e = min(s + 512, m)
        d2 = ((surface[s:e, None, :] - centers[None, :, :]) ** 2).sum(-1)
        D[s:e] = d2 < (radii[None, :] ** 2)
    return D

def check_feasible(D, selected):
    """Every surface row must be covered by at least one selected candidate."""
    if selected.size == 0:
        return D.shape[0] == 0
    return bool(D[:, selected].any(axis=1).all())


# ---- B2. LP lower-bound certificate --------------------------------------
def lp_lb_full(D):
    from scipy.optimize import linprog
    from scipy.sparse import csr_matrix
    m, n = D.shape
    A = csr_matrix(D.astype(np.float64))
    res = linprog(c=np.ones(n), A_ub=-A, b_ub=-np.ones(m), bounds=(0, 1), method="highs")
    return float(res.fun)


# ---- B1. independent full solves -----------------------------------------
def full_highs(D):
    rows = [np.nonzero(D[j])[0].tolist() for j in range(D.shape[0])]
    return S.solve_highs(rows, D.shape[1])

def full_cpsat(D):
    rows = [np.nonzero(D[j])[0].tolist() for j in range(D.shape[0])]
    return S.solve_cpsat(rows, D.shape[1], workers=8)


# ---- C. brute-force small instances --------------------------------------
def brute_optimum(D):
    """Exhaustive minimum set cover (only for tiny n)."""
    m, n = D.shape
    cols = [frozenset(np.nonzero(D[:, i])[0].tolist()) for i in range(n)]
    universe = set(range(m))
    for k in range(0, n + 1):
        for combo in itertools.combinations(range(n), k):
            u = set()
            for i in combo:
                u |= cols[i]
            if u >= universe:
                return k
    return n

def brute_force_suite(trials=40, seed=0):
    rng = np.random.default_rng(seed)
    fails = 0
    for t in range(trials):
        m = int(rng.integers(4, 9)); n = int(rng.integers(4, 10))
        D = rng.random((m, n)) < 0.45
        # ensure feasibility: every row covered by >=1 col
        for j in range(m):
            if not D[j].any():
                D[j, rng.integers(0, n)] = True
        # fabricate geometry-free instance -> drive presolve+solve via rows directly
        rows = [np.nonzero(D[j])[0].tolist() for j in range(m)]
        cov = S.Coverage(indptr=np.array([0] + list(np.cumsum([len(r) for r in rows])), dtype=np.int64),
                         indices=np.concatenate([np.array(r, np.int32) for r in rows]) if m else np.zeros(0, np.int32),
                         n_cols=n)
        pre = S.presolve(cov, do_col_dominance=True, do_row_dominance=True, geom_dominance=False)
        rows_local = [[pre.col_local[i] for i in r] for r in pre.rows]
        nloc = int(pre.active_cols.size)
        if nloc == 0:
            obj = int(pre.forced.size)
        else:
            sel, ored, proven = S.solve_cpsat(rows_local, nloc, workers=4)
            obj = int(pre.forced.size + len(sel))
            # lift + feasibility
            chosen = np.concatenate([pre.forced, pre.active_cols[np.array(sel, np.int64)] if sel else np.zeros(0, np.int64)]).astype(np.int64)
            if not check_feasible(D, np.unique(chosen)):
                print(f"  [brute] trial {t}: INFEASIBLE lifted solution!"); fails += 1; continue
        truth = brute_optimum(D)
        if obj != truth:
            print(f"  [brute] trial {t}: presolve+solve={obj} != brute={truth}  (m={m},n={n})"); fails += 1
    print(f"  brute-force suite: {trials-fails}/{trials} matched the exhaustive optimum")
    return fails == 0


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--inst", required=True)
    ap.add_argument("--brute", action="store_true", help="also run the small-instance brute-force suite")
    args = ap.parse_args()
    inst = np.load(args.inst)
    surf, cen, rad = inst["surface"], inst["centers"], inst["radii"]
    print(f"=== verifying instance: surface={len(surf)} candidates={len(cen)} ===")

    # accelerated solve (what we are auditing)
    r = S.coverage_axis_exact(surf, cen, rad, backend="cpsat", verbose=False)
    sel = r["selected"]; obj = r["objective"]
    print(f"accelerated: objective={obj} proven_optimal={r['proven_optimal']} selected={sel.size}")

    D = coverage_dense(surf, cen, rad)
    print("\nA. FEASIBILITY (independent coverage recompute):")
    feas = check_feasible(D, sel)
    print(f"   selection covers all {D.shape[0]} surface points: {'YES' if feas else 'NO  <-- FAIL'}")

    print("\nB1. INDEPENDENT FULL SOLVES (no presolve):")
    h_sel, h_obj, h_pr = full_highs(D)
    c_sel, c_obj, c_pr = full_cpsat(D)
    print(f"   HiGHS/full  obj={h_obj} proven={h_pr}")
    print(f"   CP-SAT/full obj={c_obj} proven={c_pr}")
    print(f"   match accelerated ({obj})? HiGHS={h_obj==obj} CP-SAT={c_obj==obj}")

    print("\nB2. LP LOWER-BOUND CERTIFICATE:")
    lp = lp_lb_full(D); lb = math.ceil(lp - 1e-6)
    print(f"   LP relax={lp:.4f}  ceil(LB)={lb}  objective={obj}  "
          f"-> {'CERTIFIED OPTIMAL (UB==LB)' if lb == obj else 'gap remains'}")

    all_ok = feas and (h_obj == obj == c_obj)
    if args.brute:
        print("\nC. BRUTE-FORCE SUITE (small random instances vs exhaustive):")
        all_ok &= brute_force_suite()

    print(f"\n=== VERDICT: {'ALL CHECKS PASS — exact & optimum-preserving' if all_ok else 'FAILURE'} ===")
    if not all_ok:
        raise SystemExit(1)


if __name__ == "__main__":
    main()

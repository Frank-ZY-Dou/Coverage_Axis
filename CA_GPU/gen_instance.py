"""Build a Coverage-Axis SCP instance (surface samples, inner candidates, radii)
from a watertight .off, and save it as an .npz so the baseline and the
accelerated exact solver consume byte-identical input.

Mirrors the reference pipeline: candidates = points inside the surface
(generalized winding number > 0.5), radius_i = (nearest surface distance) +
dilation. Coverage is then `||s_j - c_i|| < r_i` (built inside scp_exact)."""
from __future__ import annotations
import argparse
from pathlib import Path
import numpy as np


def parse_off(path: Path):
    """Token-stream OFF reader. Tolerates blank/comment lines and vertex rows
    with extra columns (normals/colors); reads only triangle faces."""
    toks = []
    with path.open() as f:
        for line in f:
            line = line.split("#", 1)[0].strip()
            if line:
                toks.append(line)
    # first token line is 'OFF' (optionally followed by the counts)
    first = toks[0].split()
    rest = first[1:] if first[0].upper().startswith("OFF") else first
    idx = 1
    if not rest:
        rest = toks[1].split(); idx = 2
    nv, nf = int(rest[0]), int(rest[1])
    vals = (" ".join(toks[idx:])).split()
    p = 0
    V = np.empty((nv, 3), dtype=np.float64)
    for i in range(nv):
        V[i] = [float(vals[p]), float(vals[p + 1]), float(vals[p + 2])]
        p += 3
    F = np.empty((nf, 3), dtype=np.int64)
    for i in range(nf):
        k = int(vals[p]); p += 1
        F[i] = [int(vals[p]), int(vals[p + 1]), int(vals[p + 2])]
        p += k
    return V, F


def winding_number(Q, V, F, chunk=64):
    a, b, c = V[F[:, 0]], V[F[:, 1]], V[F[:, 2]]
    out = np.empty(Q.shape[0]); inv4pi = 1.0 / (4.0 * np.pi)
    for s in range(0, Q.shape[0], chunk):
        q = Q[s:s+chunk][:, None, :]
        va, vb, vc = a[None]-q, b[None]-q, c[None]-q
        na = np.linalg.norm(va, axis=2); nb = np.linalg.norm(vb, axis=2); nc = np.linalg.norm(vc, axis=2)
        numer = (va * np.cross(vb, vc)).sum(2)
        denom = na*nb*nc + (va*vb).sum(2)*nc + (vb*vc).sum(2)*na + (vc*va).sum(2)*nb
        out[s:s+chunk] = (2.0*np.arctan2(numer, denom)).sum(1) * inv4pi
    return out


def nn_dist(P, Ref, chunk=512):
    out = np.empty(len(P))
    for s in range(0, len(P), chunk):
        d = np.sqrt(((P[s:s+chunk][:, None, :] - Ref[None, :, :]) ** 2).sum(2))
        out[s:s+chunk] = d.min(1)
    return out


def sample_surface(V, F, n, rng):
    """Area-weighted uniform samples on the triangle mesh surface."""
    area = 0.5 * np.linalg.norm(np.cross(V[F[:, 1]] - V[F[:, 0]], V[F[:, 2]] - V[F[:, 0]]), axis=1)
    fi = rng.choice(len(F), size=n, p=area / area.sum())
    u, v = rng.random(n), rng.random(n)
    over = u + v > 1.0; u[over] = 1 - u[over]; v[over] = 1 - v[over]
    a, b, c = V[F[fi, 0]], V[F[fi, 1]], V[F[fi, 2]]
    return a + u[:, None] * (b - a) + v[:, None] * (c - a)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--off", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--n-candidates", type=int, default=3000)
    ap.add_argument("--n-surface", type=int, default=3000)
    ap.add_argument("--dilation", type=float, default=0.04, help="ball dilation, fraction of bbox diagonal")
    ap.add_argument("--on-surface", action="store_true",
                    help="sample candidates ON a surface (fixed cover radius) instead of inside")
    ap.add_argument("--candidate-off", default=None,
                    help="mesh to sample candidates from (.obj/.off); defaults to --off")
    ap.add_argument("--cover-radius", type=float, default=0.08,
                    help="surface-candidate ball radius, fraction of bbox diagonal (with --on-surface)")
    ap.add_argument("--seed", type=int, default=0)
    args = ap.parse_args()
    rng = np.random.default_rng(args.seed)
    V, F = parse_off(Path(args.off))
    lo, hi = V.min(0), V.max(0); bbdiag = float(np.linalg.norm(hi - lo))
    dil = args.dilation * bbdiag

    if args.on_surface:
        # Candidates are sampled on the candidate mesh (e.g. the original
        # surface), and they must cover the points sampled on --off (e.g. the
        # double cover's two offset shells). Every ball has the same radius.
        from make_double_cover import load_mesh
        cV, cF = load_mesh(Path(args.candidate_off)) if args.candidate_off else (V, F)
        C = sample_surface(cV, cF, args.n_candidates, rng)
        S = sample_surface(V, F, args.n_surface, rng)
        R = np.full(len(C), args.cover_radius * bbdiag)
        n0 = len(S); keep = np.zeros(n0, dtype=bool)
        for s in range(0, n0, 512):
            e = min(s + 512, n0)
            d2 = ((S[s:e, None, :] - C[None, :, :]) ** 2).sum(-1)
            keep[s:e] = (d2 < (R[None, :] ** 2)).any(1)
        S = S[keep]
        np.savez(args.out, surface=S, centers=C, radii=R, bbdiag=bbdiag, dilation=R[0])
        print(f"[gen/surface] {Path(args.off).name}: surface={len(S)} (dropped {n0-len(S)}) "
              f"candidates={len(C)} cover_radius={R[0]:.4f} bbdiag={bbdiag:.4f} -> {args.out}")
        return

    cand, tries = [], 0
    while sum(len(x) for x in cand) < args.n_candidates and tries < 80:
        tries += 1
        batch = rng.uniform(lo, hi, size=(max(8000, args.n_candidates * 4), 3))
        cand.append(batch[np.abs(winding_number(batch, V, F)) > 0.5])
    C = np.concatenate(cand)[:args.n_candidates]
    R = nn_dist(C, V) + dil
    surf_idx = rng.choice(len(V), size=min(args.n_surface, len(V)), replace=False)
    S = V[surf_idx]
    # Keep only surface points coverable by at least one candidate ball
    # (||s - c_i|| < r_i for some i). Uncoverable points make the SCP infeasible;
    # the reference pipeline likewise covers only coverable samples. Reported,
    # never silent.
    n0 = len(S)
    keep = np.zeros(n0, dtype=bool)
    for s in range(0, n0, 512):
        e = min(s + 512, n0)
        d2 = ((S[s:e, None, :] - C[None, :, :]) ** 2).sum(-1)
        keep[s:e] = (d2 < (R[None, :] ** 2)).any(1)
    S = S[keep]
    np.savez(args.out, surface=S, centers=C, radii=R, bbdiag=bbdiag, dilation=dil)
    print(f"[gen] {Path(args.off).name}: surface={len(S)} (dropped {n0-len(S)} "
          f"uncoverable) candidates={len(C)} dilation={dil:.4f} bbdiag={bbdiag:.4f} -> {args.out}")


if __name__ == "__main__":
    main()

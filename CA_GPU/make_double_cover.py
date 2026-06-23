"""Build the double cover of a surface mesh by offsetting it inward and outward.

Each vertex is moved +delta and -delta along its normal, giving two parallel
shells (the orientation double cover of an orientable surface, realized as two
disjoint sheets). The two shells are written as one combined .off.

Usage:
    python make_double_cover.py --in ../input/Torus.obj --out double_torus.off --delta 0.025
"""
from __future__ import annotations
import argparse
from pathlib import Path
import numpy as np


def load_mesh(path: Path):
    """Load a triangle mesh from .obj or .off (triangulating polygon faces)."""
    V, F = [], []
    if path.suffix.lower() == ".obj":
        for ln in path.open():
            t = ln.split()
            if not t:
                continue
            if t[0] == "v":
                V.append([float(t[1]), float(t[2]), float(t[3])])
            elif t[0] == "f":
                idx = [int(p.split("/")[0]) - 1 for p in t[1:]]
                for k in range(1, len(idx) - 1):
                    F.append([idx[0], idx[k], idx[k + 1]])
        return np.array(V, np.float64), np.array(F, np.int64)
    # .off (token stream; tolerates blanks/comments)
    toks = [l.split("#", 1)[0].strip() for l in path.open()]
    toks = [t for t in toks if t]
    first = toks[0].split()
    rest = first[1:] if first[0].upper().startswith("OFF") else first
    idx = 1
    if not rest:
        rest = toks[1].split(); idx = 2
    nv, nf = int(rest[0]), int(rest[1])
    vals = (" ".join(toks[idx:])).split(); p = 0
    for _ in range(nv):
        V.append([float(vals[p]), float(vals[p + 1]), float(vals[p + 2])]); p += 3
    for _ in range(nf):
        k = int(vals[p]); p += 1
        ids = [int(vals[p + t]) for t in range(k)]; p += k
        for t in range(1, k - 1):
            F.append([ids[0], ids[t], ids[t + 1]])
    return np.array(V, np.float64), np.array(F, np.int64)


def vertex_normals(V, F):
    fn = np.cross(V[F[:, 1]] - V[F[:, 0]], V[F[:, 2]] - V[F[:, 0]])
    N = np.zeros_like(V)
    for k in range(3):
        np.add.at(N, F[:, k], fn)
    n = np.linalg.norm(N, axis=1, keepdims=True); n[n == 0] = 1.0
    return N / n


def write_off(path: Path, V, F):
    with path.open("w") as f:
        f.write("OFF\n")
        f.write(f"{len(V)} {len(F)} 0\n")
        for x, y, z in V:
            f.write(f"{x:.6f} {y:.6f} {z:.6f}\n")
        for a, b, c in F:
            f.write(f"3 {a} {b} {c}\n")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--in", dest="inp", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--delta", type=float, default=0.025, help="offset distance (mesh units)")
    args = ap.parse_args()

    V, F = load_mesh(Path(args.inp))
    N = vertex_normals(V, F)
    outer = V + args.delta * N
    inner = V - args.delta * N
    Vc = np.vstack([outer, inner])
    Fc = np.vstack([F, F[:, ::-1] + len(V)])      # inner shell wound opposite
    write_off(Path(args.out), Vc, Fc)
    print(f"[double-cover] {Path(args.inp).name}: {len(V)} verts -> {len(Vc)} "
          f"(two shells, +/-{args.delta}) -> {args.out}")


if __name__ == "__main__":
    main()

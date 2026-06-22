"""Polyscope visualization of the Coverage-Axis SCP selection.

Solves the same instance two ways and renders each over the input surface, with
the selected inner points (Coverage Axis poles) drawn as uniform small spheres:

  * reference : HiGHS milp on the full coverage matrix (as in Coverage_Axis_mesh.py)
  * gpu       : the GPU-accelerated exact solver (scp_fast)

Both are exact, so the objective is identical; the figures let you see that the
GPU version selects an equally-optimal pole set. Headless-friendly (EGL); writes
`<name>_reference_scp.png` and `<name>_gpu_scp.png` to figs/.

Usage:
    python gen_instance.py --off ../input/01Ants-12_mesh.off --out ants.npz \
        --n-candidates 3000 --n-surface 3000 --dilation 0.03
    python viz_scp.py --off ../input/01Ants-12_mesh.off --inst ants.npz --name ants
"""
from __future__ import annotations

import argparse
from pathlib import Path
import numpy as np
import polyscope as ps

import scp_exact as SX
import scp_fast as SF
from gen_instance import parse_off


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--off", required=True, help="input surface mesh (.off)")
    ap.add_argument("--inst", required=True, help="SCP instance .npz from gen_instance.py")
    ap.add_argument("--name", default="mesh")
    ap.add_argument("--out", default="figs")
    args = ap.parse_args()
    out_dir = Path(args.out).resolve(); out_dir.mkdir(parents=True, exist_ok=True)

    surf_v, surf_f = parse_off(Path(args.off))
    inst = np.load(args.inst)
    surface, centers, radii = inst["surface"], inst["centers"], inst["radii"]

    # reference solver: HiGHS milp on the full coverage matrix
    cov = SX.build_coverage(surface, centers, radii)
    rows = [list(int(i) for i in cov.row(j)) for j in range(cov.m)]
    selR, objR, prR = SX.solve_highs(rows, cov.n_cols)
    selR = np.array(selR, dtype=np.int64)
    # GPU-accelerated exact solver
    resG = SF.coverage_axis_fast(surface, centers, radii, backend="cpsat", verbose=False)
    selG = resG["selected"]
    identical = set(selR.tolist()) == set(selG.tolist())
    print(f"{args.name}: reference obj={objR} (proven={prR}) | gpu obj={resG['objective']} "
          f"(proven={resG['proven_optimal']}) | same selection={identical}"
          f"{'' if identical else ' (equally-optimal, different pole set)'}")

    ps.set_allow_headless_backends(True)
    ps.set_program_name("Coverage Axis SCP viz")
    ps.set_window_size(1280, 960)
    ps.set_ground_plane_mode("none")
    ps.set_SSAA_factor(2)
    ps.init()

    bb_min, bb_max = surf_v.min(0), surf_v.max(0)
    diag = float(np.linalg.norm(bb_max - bb_min)); ctr = (bb_min + bb_max) * 0.5
    uniform = 0.008 * diag           # selected-pole sphere radius
    # Side profile from slightly above: look along the lateral (mid) bbox axis so
    # the longest axis (body) is horizontal, "up" the shortest axis, camera raised.
    e = np.eye(3)
    order = np.argsort(-(bb_max - bb_min))
    lateral_a, up_a = e[order[1]], e[order[2]]
    cam = ctr + lateral_a * 1.7 * diag + up_a * 0.85 * diag    # 15% closer

    ps.register_surface_mesh("input surface", surf_v, surf_f, smooth_shade=True,
                             color=(0.78, 0.78, 0.80), edge_width=0.0, transparency=0.55)

    for tag, sel, col in [("reference", selR, (0.90, 0.20, 0.20)),   # original: red
                          ("gpu", selG, (0.15, 0.45, 0.95))]:        # GPU: blue
        cloud = ps.register_point_cloud("selected Coverage-Axis poles", centers[sel],
                                        color=col, radius=uniform)
        ps.look_at_dir(camera_location=cam, target=ctr, up_dir=up_a)
        ps.frame_tick()
        out = out_dir / f"{args.name}_{tag}_scp.png"
        ps.screenshot(str(out), transparent_bg=False)
        print(f"  saved {out}  ({sel.size} poles)")
        cloud.remove()


if __name__ == "__main__":
    main()

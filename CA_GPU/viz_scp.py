"""Polyscope visualization of the Coverage-Axis SCP selection.

Two layouts (headless-friendly EGL; PNGs to figs/):

* default: render the input surface twice, with the reference (red) and GPU
  (blue) selected poles, to show the two exact solvers agree.
      python viz_scp.py --off ants.off --inst ants.npz --name ants

* with --double-off: a two-panel input/result figure for the double-cover
  example. Left = original mesh (yellow) + the double cover; right = original
  mesh (yellow) + the selected points.
      python viz_scp.py --off ../input/Torus.obj --double-off double_torus.off \
          --inst double_torus.npz --name double_torus
"""
from __future__ import annotations

import argparse
from pathlib import Path
import numpy as np
import polyscope as ps

import scp_exact as SX
import scp_fast as SF
from make_double_cover import load_mesh

YELLOW = (0.95, 0.78, 0.15)
GREY = (0.78, 0.78, 0.80)
RED = (0.90, 0.20, 0.20)
BLUE = (0.15, 0.45, 0.95)


def setup_camera(bb_min, bb_max):
    diag = float(np.linalg.norm(bb_max - bb_min)); ctr = (bb_min + bb_max) * 0.5
    e = np.eye(3)
    order = np.argsort(-(bb_max - bb_min))          # [long, mid(lateral), short(up)]
    lateral_a, up_a = e[order[1]], e[order[2]]
    cam = ctr + lateral_a * 1.7 * diag + up_a * 0.85 * diag
    return cam, ctr, up_a, diag


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--off", required=True, help="original surface mesh (.obj/.off)")
    ap.add_argument("--inst", required=True, help="SCP instance .npz from gen_instance.py")
    ap.add_argument("--double-off", default=None, help="double-cover mesh for the input panel")
    ap.add_argument("--name", default="mesh")
    ap.add_argument("--out", default="figs")
    args = ap.parse_args()
    out_dir = Path(args.out).resolve(); out_dir.mkdir(parents=True, exist_ok=True)

    surf_v, surf_f = load_mesh(Path(args.off))
    inst = np.load(args.inst)
    surface, centers, radii = inst["surface"], inst["centers"], inst["radii"]

    cov = SX.build_coverage(surface, centers, radii)
    rows = [list(int(i) for i in cov.row(j)) for j in range(cov.m)]
    selR, objR, prR = SX.solve_highs(rows, cov.n_cols)
    selR = np.array(selR, dtype=np.int64)
    resG = SF.coverage_axis_fast(surface, centers, radii, backend="cpsat", verbose=False)
    selG = resG["selected"]
    print(f"{args.name}: reference obj={objR} (proven={prR}) | gpu obj={resG['objective']} "
          f"(proven={resG['proven_optimal']})")

    ps.set_allow_headless_backends(True)
    ps.set_program_name("Coverage Axis viz")
    ps.set_window_size(1280, 960)
    ps.set_ground_plane_mode("none")
    ps.set_SSAA_factor(2)
    ps.init()

    cam, ctr, up_a, diag = setup_camera(surf_v.min(0), surf_v.max(0))
    uniform = 0.008 * diag

    def shot(tag):
        ps.look_at_dir(camera_location=cam, target=ctr, up_dir=up_a)
        ps.frame_tick()
        out = out_dir / f"{args.name}_{tag}.png"
        ps.screenshot(str(out), transparent_bg=False)
        print(f"  saved {out}")

    if args.double_off:
        dv, df = load_mesh(Path(args.double_off))
        # left: original (yellow) + double cover (grey, translucent)
        om = ps.register_surface_mesh("original surface", surf_v, surf_f, smooth_shade=True,
                                      color=YELLOW, edge_width=0.0, transparency=0.85)
        dc = ps.register_surface_mesh("double cover", dv, df, smooth_shade=True,
                                      color=GREY, edge_width=0.0, transparency=0.35)
        shot("input")
        # right: original (yellow) + selected points (blue)
        dc.set_enabled(False)
        om.set_transparency(0.6)
        pts = ps.register_point_cloud("selected points", centers[selG], color=BLUE, radius=uniform)
        shot("result")
        pts.remove()
    else:
        ps.register_surface_mesh("input surface", surf_v, surf_f, smooth_shade=True,
                                 color=GREY, edge_width=0.0, transparency=0.55)
        for tag, sel, col in [("reference_scp", selR, RED), ("gpu_scp", selG, BLUE)]:
            cloud = ps.register_point_cloud("selected Coverage-Axis poles", centers[sel],
                                            color=col, radius=uniform)
            shot(tag)
            cloud.remove()


if __name__ == "__main__":
    main()

# Author: Frank ZY Dou
"""Skeleton connection for a mesh input with Q-MAT (Sec. 4.3 of the Coverage Axis paper, "mesh with candidates from
Voronoi diagram").

The initial medial axis of the mesh (the inner Voronoi diagram written by mesh_to_ma.py) is simplified by the edge
collapse of Q-MAT while the inner points selected by Coverage Axis / Coverage Axis++ serve as anchors: an edge whose
two endpoints are anchors is never collapsed, an edge with one anchor endpoint is collapsed onto the anchor, and an
edge without anchor endpoints is collapsed onto the optimal sphere of Q-MAT, until only the anchors remain. This runs
mode 2 of the Q-MAT executable built from skel_connection/QMAT (every selected point is first snapped to its nearest
medial axis vertex, an exact match for the Voronoi candidates).

Q-MAT works in coordinates divided by the bounding-box diagonal of the mesh and compares the selected points with
its vertices in that scale, so the points are written for it in that scale; the exported skeleton is in the original
coordinates again.

Outputs ``<out>.obj`` (vertices ``v``, edges ``l``, triangles ``f``), ``<out>.ma`` (Q-MAT format with radii) and
``<out>_balls.obj`` (the medial balls of the skeleton vertices as sphere meshes, for visualization).

Example:
    python Coverage_Axis_connection.py --mesh ./input/01Ants-12_mesh.off --ma ./input/01Ants-12_mesh.ma \
        --selected ./output/mesh_selected_inner_points.txt --out ./output/mesh_skeleton
"""
import argparse
import glob
import os
import shutil
import subprocess
import tempfile
import numpy as np
import trimesh
from scipy.spatial import cKDTree
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components
from utils import read_VD


def read_ma(path):
    with open(path) as f:
        f.readline()
        verts, edges, faces = [], [], []
        for line in f:
            parts = line.split()
            if not parts:
                continue
            if parts[0] == 'v':
                verts.append([float(x) for x in parts[1:5]])
            elif parts[0] == 'e':
                edges.append([int(parts[1]), int(parts[2])])
            elif parts[0] == 'f':
                faces.append([int(parts[1]), int(parts[2]), int(parts[3])])
    return np.array(verts), np.array(edges, dtype=int).reshape(-1, 2), np.array(faces, dtype=int).reshape(-1, 3)


def write_obj(path, verts, edges, faces):
    with open(path, 'w') as f:
        for v in verts:
            f.write('v %f %f %f\n' % (v[0], v[1], v[2]))
        for e in edges:
            f.write('l %d %d\n' % (e[0] + 1, e[1] + 1))
        for t in faces:
            f.write('f %d %d %d\n' % (t[0] + 1, t[1] + 1, t[2] + 1))


def write_balls(path, verts, sphere='./assets/sphere_I.obj'):
    """The medial balls of the skeleton vertices as one OBJ of sphere meshes."""
    unit = trimesh.load(sphere)
    base_verts, base_faces = np.asarray(unit.vertices), np.asarray(unit.faces)
    with open(path, 'w') as f:
        for v in verts:
            for q in base_verts * v[3] + v[:3]:
                f.write('v %f %f %f\n' % (q[0], q[1], q[2]))
        for i in range(len(verts)):
            for t in base_faces + i * len(base_verts) + 1:
                f.write('f %d %d %d\n' % (t[0], t[1], t[2]))


def main():
    parser = argparse.ArgumentParser(description='Connect the selected inner points of a mesh into a skeleton with Q-MAT.')
    parser.add_argument('--mesh', required=True, help='input triangle mesh (.off)')
    parser.add_argument('--ma', required=True, help='initial medial axis of the mesh (.ma, see mesh_to_ma.py)')
    parser.add_argument('--selected', required=True, help='selected inner points with radii (v x y z r per line)')
    parser.add_argument('--qmat', default='./skel_connection/QMAT/build/QMAT', help='Q-MAT executable (default: ./skel_connection/QMAT/build/QMAT)')
    parser.add_argument('--out', required=True, help='output prefix; writes <out>.obj and <out>.ma')
    args = parser.parse_args()

    if not os.path.isfile(args.qmat):
        raise FileNotFoundError("Q-MAT executable not found at %s; build it as described in the README (Installation)" % args.qmat)
    mesh = trimesh.load(args.mesh)
    diagonal = float(np.linalg.norm(mesh.bounds[1] - mesh.bounds[0]))
    selected, radius = read_VD(args.selected)
    selected = np.concatenate([np.array(selected, dtype=float), np.array(radius, dtype=float)], axis=1)
    print("Selected inner points: %d" % len(selected))

    workdir = tempfile.mkdtemp(prefix='qmat_')
    poles = os.path.join(workdir, 'selected_points.txt')
    with open(poles, 'w') as f:
        for p in selected / diagonal:
            f.write('v %.12f %.12f %.12f %.12f\n' % (p[0], p[1], p[2], p[3]))
    command = [os.path.abspath(args.qmat), '2', os.path.abspath(args.mesh), os.path.abspath(args.ma), '0', workdir + '/', poles]
    print("Running: " + ' '.join(command))
    result = subprocess.run(command, stdin=subprocess.DEVNULL, capture_output=True, text=True)
    exported = glob.glob(os.path.join(workdir, 'export_half___v_*.ma'))
    if result.returncode != 0 or not exported:
        raise RuntimeError("Q-MAT failed (exit code %d):\n%s\n%s" % (result.returncode, result.stdout[-2000:], result.stderr[-2000:]))
    verts, edges, faces = read_ma(exported[0])
    shutil.copyfile(exported[0], args.out + '.ma')
    shutil.rmtree(workdir)

    n = len(verts)
    adjacency = coo_matrix((np.ones(len(edges)), (edges[:, 0], edges[:, 1])), shape=(n, n))
    n_components = connected_components(adjacency, directed=False)[0]
    distance = cKDTree(verts[:, :3]).query(selected[:, :3])[0]
    print("Skeleton: %d vertices, %d edges, %d triangles, %d connected components" % (n, len(edges), len(faces), n_components))
    print("Distance from the selected points to the skeleton vertices: mean %.4f, max %.4f (bounding-box diagonal %.3f)"
          % (distance.mean(), distance.max(), diagonal))

    write_obj(args.out + '.obj', verts, edges, faces)
    write_balls(args.out + '_balls.obj', verts, os.path.join(os.path.dirname(os.path.abspath(__file__)), 'assets', 'sphere_I.obj'))
    print("Saved %s.obj, %s.ma and %s_balls.obj" % (args.out, args.out, args.out))


if __name__ == '__main__':
    main()

# Author: Frank ZY Dou
"""Initial medial axis of a mesh for Q-MAT: the Voronoi diagram of the mesh vertices restricted to the inside of the mesh.

The vertices of the medial mesh are the circumcentres of the Delaunay tetrahedra of the mesh vertices that lie inside the
mesh (generalized winding number > 0.5), with the circumradius as ball radius. Two vertices are connected when their
tetrahedra share a triangle (a Voronoi edge), and the Voronoi face dual to a Delaunay edge whose surrounding tetrahedra are
all inside is written as a fan of triangles. The result is written in the .ma format of Q-MAT (0-based indices):

    nv ne nf
    v x y z r
    e i j
    f i j k

Example:
    python mesh_to_ma.py --mesh ./input/01Ants-12_mesh.off --out ./input/01Ants-12_mesh.ma
"""
import argparse
import numpy as np
import torch
import trimesh
from scipy.spatial import Delaunay
from utils import winding_number


def circumspheres(points, tets):
    a = points[tets[:, 0]]
    m = np.stack([points[tets[:, 1]] - a, points[tets[:, 2]] - a, points[tets[:, 3]] - a], axis=1)
    rhs = 0.5 * (m ** 2).sum(-1)
    x = np.linalg.solve(m, rhs[..., None])[..., 0]
    return a + x, np.linalg.norm(x, axis=1)


def inside(points, vertices, faces, batch=5000):
    v = torch.tensor(vertices).cuda().double()
    f = torch.tensor(faces).cuda().long()
    w = [winding_number(torch.tensor(points[i:i + batch]).cuda().double(), v, f) for i in range(0, len(points), batch)]
    return torch.cat(w).cpu().numpy() > 0.5


def inner_voronoi_diagram(points, faces):
    dt = Delaunay(points)
    tets, neighbors = dt.simplices, dt.neighbors
    centers, radii = circumspheres(points, tets)
    keep = inside(centers, points, faces)
    index = np.full(len(tets), -1)
    index[keep] = np.arange(keep.sum())

    # Voronoi edges: dual to the triangles shared by two inner tetrahedra
    t = np.repeat(np.arange(len(tets)), 4)
    n = neighbors.reshape(-1)
    ok = (n > t) & keep[t] & keep[np.maximum(n, 0)]
    edges = np.stack([index[t[ok]], index[n[ok]]], axis=1)

    # Voronoi faces: dual to the Delaunay edges that are not on the convex hull and whose tetrahedra are all inner
    pairs = np.array([[0, 1], [0, 2], [0, 3], [1, 2], [1, 3], [2, 3]])
    te = np.sort(tets[:, pairs], axis=2).reshape(-1, 2)           # (n_tets * 6, 2): the edges of each tetrahedron
    te_tet = np.repeat(np.arange(len(tets)), 6)
    dedges, inv = np.unique(te, axis=0, return_inverse=True)
    inv = inv.reshape(-1)
    hull = np.sort(dt.convex_hull[:, [[0, 1], [0, 2], [1, 2]]], axis=2).reshape(-1, 2)
    on_hull = np.isin(dedges[:, 0] * len(points) + dedges[:, 1], hull[:, 0] * len(points) + hull[:, 1])
    all_inner = np.ones(len(dedges), dtype=bool)
    np.logical_and.at(all_inner, inv, keep[te_tet])
    good = all_inner & ~on_hull
    order = np.argsort(inv, kind='stable')
    inv, te_tet = inv[order], te_tet[order]
    faces, diagonals = [], []
    starts = np.r_[0, np.flatnonzero(np.diff(inv)) + 1, len(inv)]
    for s, e in zip(starts[:-1], starts[1:]):
        d = inv[s]
        if not good[d]:
            continue
        ring = te_tet[s:e]
        c = centers[ring]
        axis = points[dedges[d, 1]] - points[dedges[d, 0]]
        axis /= np.linalg.norm(axis)
        u = np.cross(axis, [1.0, 0.0, 0.0] if abs(axis[0]) < 0.9 else [0.0, 1.0, 0.0])
        u /= np.linalg.norm(u)
        w = np.cross(axis, u)
        rel = c - c.mean(0)
        ring = ring[np.argsort(np.arctan2(rel @ w, rel @ u))]
        poly = index[ring]
        for k in range(1, len(poly) - 1):
            faces.append([poly[0], poly[k], poly[k + 1]])
            if k > 1:
                diagonals.append([poly[0], poly[k]])
    faces = np.array(faces, dtype=int).reshape(-1, 3)
    if diagonals:
        edges = np.vstack([edges, np.array(diagonals, dtype=int)])
    edges = np.unique(np.sort(edges, axis=1), axis=0)
    return centers[keep], radii[keep], edges, faces


def write_ma(path, centers, radii, edges, faces):
    with open(path, 'w') as f:
        f.write('%d %d %d\n' % (len(centers), len(edges), len(faces)))
        for c, r in zip(centers, radii):
            f.write('v %.15f %.15f %.15f %.15f\n' % (c[0], c[1], c[2], r))
        for e in edges:
            f.write('e %d %d\n' % (e[0], e[1]))
        for t in faces:
            f.write('f %d %d %d\n' % (t[0], t[1], t[2]))


def main():
    parser = argparse.ArgumentParser(description='Initial medial axis (inner Voronoi diagram of the mesh vertices) in the .ma format of Q-MAT.')
    parser.add_argument('--mesh', required=True, help='input triangle mesh (.off/.obj)')
    parser.add_argument('--out', required=True, help='output .ma file')
    args = parser.parse_args()
    mesh = trimesh.load(args.mesh)
    points, faces = np.asarray(mesh.vertices, dtype=float), np.asarray(mesh.faces)
    centers, radii, edges, tris = inner_voronoi_diagram(points, faces)
    write_ma(args.out, centers, radii, edges, tris)
    print("Mesh: %d vertices, %d faces" % (len(points), len(faces)))
    print("Medial axis: %d vertices, %d edges, %d faces -> %s" % (len(centers), len(edges), len(tris), args.out))


if __name__ == '__main__':
    main()

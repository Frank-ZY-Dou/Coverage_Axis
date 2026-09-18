# Author: Frank ZY Dou
"""Skeleton connection for a mesh input (Sec. 4.3 of the Coverage Axis paper, mesh with candidates from the Voronoi
diagram).

The initial medial axis of the mesh (its inner Voronoi diagram, written by mesh_to_ma.py) carries the correct
connectivity of the shape, and the selected inner points are vertices of it. The medial axis is simplified onto the
selected points: every vertex of the medial axis is merged into the selected point that is closest to it along the
medial axis (geodesic distance on the medial mesh, computed by a multi-source Dijkstra search). Two selected points
are connected by an edge when their regions share an edge of the medial axis, and three selected points form a
triangle when their regions meet at a face of the medial axis. The skeleton thus inherits the connectivity of the
medial axis and has the selected points as its only vertices.

Outputs ``<out>.obj`` (vertices ``v``, edges ``l``, triangles ``f``) and ``<out>.ma`` (the .ma format of Q-MAT:
``nv ne nf`` header, ``v x y z r``, 0-based ``e i j`` and ``f i j k``).

Example:
    python Coverage_Axis_connection.py --ma ./input/01Ants-12_mesh.ma --selected ./output/mesh_selected_inner_points.txt --out ./output/mesh_skeleton
"""
import argparse
import numpy as np
from scipy.spatial import cKDTree
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components, dijkstra
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


def write_ma(path, verts, edges, faces):
    with open(path, 'w') as f:
        f.write('%d %d %d\n' % (len(verts), len(edges), len(faces)))
        for v in verts:
            f.write('v %.8f %.8f %.8f %.8f\n' % (v[0], v[1], v[2], v[3]))
        for e in edges:
            f.write('e %d %d\n' % (e[0], e[1]))
        for t in faces:
            f.write('f %d %d %d\n' % (t[0], t[1], t[2]))


def connect(ma_verts, ma_edges, ma_faces, anchors):
    """Merge every medial axis vertex into its nearest anchor along the medial axis and return the edges and
    triangles between the regions (indices into `anchors`), and the number of medial axis vertices no anchor reaches."""
    n = len(ma_verts)
    length = np.linalg.norm(ma_verts[ma_edges[:, 0], :3] - ma_verts[ma_edges[:, 1], :3], axis=1) + 1e-12
    graph = coo_matrix((np.concatenate([length, length]),
                        (np.concatenate([ma_edges[:, 0], ma_edges[:, 1]]), np.concatenate([ma_edges[:, 1], ma_edges[:, 0]]))),
                       shape=(n, n)).tocsr()
    distance, _, source = dijkstra(graph, directed=False, indices=anchors, min_only=True, return_predecessors=True)
    reached = np.isfinite(distance)
    region = np.full(n, -1)
    region[anchors] = np.arange(len(anchors))
    label = np.full(n, -1)
    label[reached] = region[source[reached]]

    edge_labels = label[ma_edges]
    edge_labels = edge_labels[(edge_labels >= 0).all(1) & (edge_labels[:, 0] != edge_labels[:, 1])]
    edges = np.unique(np.sort(edge_labels, axis=1), axis=0)
    face_labels = label[ma_faces]
    face_labels = face_labels[(face_labels >= 0).all(1)]
    distinct = (face_labels[:, 0] != face_labels[:, 1]) & (face_labels[:, 1] != face_labels[:, 2]) & (face_labels[:, 0] != face_labels[:, 2])
    faces = np.unique(np.sort(face_labels[distinct], axis=1), axis=0)
    return edges, faces, int(np.sum(~reached))


def main():
    parser = argparse.ArgumentParser(description='Connect the selected inner points of a mesh into a skeleton along its medial axis.')
    parser.add_argument('--ma', required=True, help='initial medial axis of the mesh (.ma written by mesh_to_ma.py)')
    parser.add_argument('--selected', required=True, help='selected inner points with radii (v x y z r per line)')
    parser.add_argument('--out', required=True, help='output prefix; writes <out>.obj and <out>.ma')
    args = parser.parse_args()

    ma_verts, ma_edges, ma_faces = read_ma(args.ma)
    selected, radius = read_VD(args.selected)
    selected = np.concatenate([np.array(selected, dtype=float), np.array(radius, dtype=float)], axis=1)
    print("Medial axis: %d vertices, %d edges, %d faces" % (len(ma_verts), len(ma_edges), len(ma_faces)))
    print("Selected inner points: %d" % len(selected))

    snap, anchors = cKDTree(ma_verts[:, :3]).query(selected[:, :3])
    print("Distance from the selected points to their medial axis vertices: max %.6f" % snap.max())
    if len(np.unique(anchors)) != len(anchors):
        raise ValueError("Two selected points fall on the same medial axis vertex. Select among the vertices of the medial axis "
                         "(inner_points = \"voronoi\" with the candidates written by mesh_to_ma.py).")

    edges, faces, unreached = connect(ma_verts, ma_edges, ma_faces, anchors)
    if unreached:
        print("Warning: %d medial axis vertices belong to parts of the medial axis without any selected point" % unreached)

    n = len(selected)
    adjacency = coo_matrix((np.ones(len(edges)), (edges[:, 0], edges[:, 1])), shape=(n, n))
    n_components = connected_components(adjacency, directed=False)[0]
    degree = np.bincount(edges.ravel(), minlength=n)
    print("Skeleton: %d vertices, %d edges, %d triangles, %d connected components, %d isolated vertices"
          % (n, len(edges), len(faces), n_components, int(np.sum(degree == 0))))

    write_obj(args.out + '.obj', selected, edges, faces)
    write_ma(args.out + '.ma', selected, edges, faces)
    print("Saved %s.obj and %s.ma" % (args.out, args.out))


if __name__ == '__main__':
    main()

import trimesh
def read_VD(path):
    points = []
    radius = []
    with open(path, "r") as f:
        for line in f.readlines():
            line = line.strip('\n')  # 去掉列表中每一个元素的换行符
            line = line.split(' ')
            points.append([float(line[1]), float(line[2]), float(line[3])])
            radius.append([float(line[4])])
    return points, radius


def save_obj(path, verts, faces=None):
    verts = verts.tolist()

    with open(path, 'w') as f:
        for v in verts:
            f.write('v %f %f %f\n' %(v[0], v[1], v[2]))
        if faces is not None:
            faces = faces.tolist()
            for ff in faces:
                f.write('f %d %d %d\n' % (ff[0] + 1, ff[1] + 1, ff[2] + 1))


# Codes borrowed from https://gist.github.com/dendenxu/ee5008acb5607195582e7983a384e644#file-moller_trumbore-py-L27
import torch



def multi_indexing(index: torch.Tensor, shape: torch.Size, dim=-2):
    shape = list(shape)
    back_pad = len(shape) - index.ndim
    for _ in range(back_pad):
        index = index.unsqueeze(-1)
    expand_shape = shape
    expand_shape[dim] = -1
    return index.expand(*expand_shape)

def multi_gather(values: torch.Tensor, index: torch.Tensor, dim=-2):
    # take care of batch dimension of, and acts like a linear indexing in the target dimention
    # we assume that the index's last dimension is the dimension to be indexed on
    return values.gather(dim, multi_indexing(index, values.shape, dim))

def winding_number(pts: torch.Tensor, verts: torch.Tensor, faces: torch.Tensor) -> torch.Tensor:
    """
    Parallel implementation of the Generalized Winding Number of points on the mesh
    O(n_points * n_faces) memory usage, parallelized execution
    1. Project tris onto the unit sphere around every points
    2. Compute the signed solid angle of the each triangle for each point
    3. Sum the solid angle of each triangle
    Parameters
    ----------
    pts    : torch.Tensor, (n_points, 3)
    verts  : torch.Tensor, (n_verts, 3)
    faces  : torch.Tensor, (n_faces, 3)
    This implementation is also able to take a/multiple batch dimension
    """
    # projection onto unit sphere: verts implementation gives a little bit more performance
    uv = verts[..., None, :, :] - pts[..., :, None, :]  # n_points, n_verts, 3
    uv = uv / uv.norm(dim=-1, keepdim=True)  # n_points, n_verts, 3

    # gather from the computed vertices (will result in a copy for sure)
    expanded_faces = faces[..., None, :, :].expand(*faces.shape[:-2], pts.shape[-2], *faces.shape[-2:])  # n_points, n_faces, 3

    u0 = multi_gather(uv, expanded_faces[..., 0])  # n, f, 3
    u1 = multi_gather(uv, expanded_faces[..., 1])  # n, f, 3
    u2 = multi_gather(uv, expanded_faces[..., 2])  # n, f, 3

    e0 = u1 - u0  # n, f, 3
    e1 = u2 - u1  # n, f, 3
    del u1

    # compute solid angle signs
    sign = (torch.cross(e0, e1) * u2).sum(dim=-1).sign()

    e2 = u0 - u2
    del u0, u2

    l0 = e0.norm(dim=-1)
    del e0

    l1 = e1.norm(dim=-1)
    del e1

    l2 = e2.norm(dim=-1)
    del e2

    # compute edge lengths: pure triangle
    l = torch.stack([l0, l1, l2], dim=-1)  # n_points, n_faces, 3

    # compute spherical edge lengths
    l = 2 * (l/2).arcsin()  # n_points, n_faces, 3

    # compute solid angle: preparing: n_points, n_faces
    s = l.sum(dim=-1) / 2
    s0 = s - l[..., 0]
    s1 = s - l[..., 1]
    s2 = s - l[..., 2]

    # compute solid angle: and generalized winding number: n_points, n_faces
    eps = 1e-10  # NOTE: will cause nan if not bigger than 1e-10
    solid = 4 * (((s/2).tan() * (s0/2).tan() * (s1/2).tan() * (s2/2).tan()).abs() + eps).sqrt().arctan()
    signed_solid = solid * sign  # n_points, n_faces

    winding = signed_solid.sum(dim=-1) / (4 * torch.pi)  # n_points

    return winding


def ball_vis(position, r):
    for i in range(r.shape[0]):
        mesh = trimesh.load('./assets/sphere_I.obj')
        mesh.vertices = mesh.vertices * r[i]
        mesh.vertices = mesh.vertices + position[i]
        mesh.export('./vis_ball/%04d.obj'%i)




# Author: Frank ZY Dou

import numpy as np
import trimesh
from scipy.optimize import milp, Bounds, LinearConstraint
from utils import  save_obj,read_VD
import torch


real_name = '01Ants-12'
surface_sample_num = 2000
dilation = 0.02
max_time_SCP = 90
# inner_points = "voronoi"
inner_points = "random"

mesh = trimesh.load('./input/%s.off'%real_name)
point_set = trimesh.sample.sample_surface(mesh, surface_sample_num)

if inner_points == "voronoi":
    medial_path ='_VD.txt'
    inner_point_path  = './input/'+real_name+medial_path
    inner_points, radius = read_VD(inner_point_path)
    inner_points = np.array(inner_points)
    radius = np.array(radius)
    radius = radius + dilation

else:
    print("Generating random samples inside the shape...")
    inner_points = trimesh.sample.volume_mesh(mesh, 30000)
    '''
    Note the returned points can be less than requested 10000 samples
    '''
    assert len(inner_points) > surface_sample_num * 4, "The number of samples: %d at least 4 times of surface samples"%len(inner_points)
    inner_points = np.array(inner_points)
    print(inner_points.shape)
    input(22)

mesh_faces = np.array(mesh.faces)
mesh_vertices = np.array(mesh.vertices)
point_set = np.array(point_set[0])

save_obj("./output/mesh.obj", mesh_vertices, mesh_faces)
save_obj("./output/mesh_samples_%d.obj"%surface_sample_num, point_set)
save_obj("./output/inner_points.obj", inner_points)

# Coverage Matrix -> GPU.
mesh_vertices_g = torch.tensor(point_set).cuda().double()
innerpoints_g = torch.tensor(inner_points).cuda().double()
radius_g = torch.tensor(radius).cuda().double()
radius_g = radius_g[:,0]
radius_g = radius_g.unsqueeze(0).repeat(len(point_set), 1)
D = torch.cdist(mesh_vertices_g, innerpoints_g, p=2)
D = torch.gt(radius_g, D).type(torch.int)
D = D.cpu().numpy()
c = np.ones(len(inner_points))
options = {"disp": True, "time_limit": max_time_SCP, }
A,b =  D, np.ones(len(point_set))
integrality = np.ones(len(inner_points))
lb, ub = np.zeros(len(inner_points)), np.ones(len(inner_points))
variable_bounds = Bounds(lb, ub)
constraints = LinearConstraint(A, lb=b)
res_milp = milp(
    c,
    integrality=integrality,
    bounds=variable_bounds,
    constraints=constraints,
    options=options)

res_milp.x = [int(x_i) for x_i in res_milp.x]
print(res_milp)
print(np.sum(res_milp.x))
value_pos = np.nonzero(res_milp.x)[0]
print("The number of selected inner points: ", len(value_pos))
save_obj("./output/selected_inner_points.obj", inner_points[value_pos])





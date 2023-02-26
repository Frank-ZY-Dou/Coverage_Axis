# Author: Frank ZY Dou
import os
import torch
import trimesh
import numpy as np
from tqdm import tqdm
from utils import  save_obj,read_VD, winding_number
from scipy.optimize import milp, Bounds, LinearConstraint


real_name = '01Ants-12_mesh'
surface_sample_num = 2000
dilation = 0.025
# inner_points = "voronoi"
inner_points = "random"
max_time_SCP = 1000 # in second


mesh = trimesh.load('./input/%s.off'%real_name)
point_set = trimesh.sample.sample_surface(mesh, surface_sample_num)

mesh_faces = np.array(mesh.faces)
mesh_vertices = np.array(mesh.vertices)
point_set = np.array(point_set[0])

if inner_points == "voronoi":
    medial_path ='_VD.txt'
    inner_point_path  = './input/'+real_name+medial_path
    inner_points, radius = read_VD(inner_point_path)
    inner_points = np.array(inner_points)
    radius = np.array(radius)
    radius = radius + dilation

else:
    print("Generating random samples inside the shape...")
    if os.path.exists("./input/%s_random.obj"%real_name):
        inner_points = trimesh.load("./input/%s_random.obj"%real_name)
        inner_points = np.array(inner_points.vertices)
        print("The number of sampled inner candidates: ", len(inner_points))
    else:
        print("Randomly Generating inner candidates...")
        random_sample_number = 500000
        min_corner = np.amin(np.amin(mesh_vertices, axis=0), axis=0)
        max_corner = np.amax(np.amax(mesh_vertices, axis=0), axis=0)
        center = min_corner + max_corner /2
        min_x = np.min(mesh_vertices[:,0])
        min_y = np.min(mesh_vertices[:,1])
        min_z = np.min(mesh_vertices[:,2])
        max_x = np.max(mesh_vertices[:, 0])
        max_y = np.max(mesh_vertices[:, 1])
        max_z = np.max(mesh_vertices[:, 2])
        P_x = (max_x) * np.random.random((random_sample_number, 1)) * 1.3 + min_x - 0.1
        P_y = (max_y) * np.random.random((random_sample_number, 1)) *1.3+ min_y - 0.1
        P_z = (max_z) * np.random.random((random_sample_number, 1))*1.3 + min_z - 0.1
        P = np.concatenate((P_x, P_y, P_z), axis=1)
        winding_con = []
        for i in tqdm(range(0, len(P), 5000)):
            start = i
            end = i + 5000
            winding = winding_number(torch.tensor(P[start:end,:]).cuda().double(), torch.tensor(mesh_vertices).cuda().double(), torch.tensor(mesh_faces).cuda().long())
            winding_con.append(winding)
        winding_con = torch.cat(winding_con, dim=0)
        inner_points = P[winding_con.cpu().numpy() > 0.5]
        save_obj("./input/%s_random.obj"%real_name, inner_points)

    inner_points_g = torch.tensor(inner_points).cuda().double()
    point_set_g = torch.tensor(point_set).cuda().double()
    dist = torch.cdist(inner_points_g, point_set_g, p=2)
    radius = dist.topk(1, largest=False).values
    radius = radius + dilation




save_obj("./output/mesh.obj", mesh_vertices, mesh_faces)
save_obj("./output/mesh_samples_%d.obj"%surface_sample_num, point_set)
save_obj("./output/mesh_inner_points.obj", inner_points)

# Coverage Matrix -> GPU.
point_set_g = torch.tensor(point_set).cuda().double()
innerpoints_g = torch.tensor(inner_points).cuda().double()
radius_g = torch.tensor(radius).cuda().double()
radius_g = radius_g[:,0]
radius_g = radius_g.unsqueeze(0).repeat(len(point_set), 1)
D = torch.cdist(point_set_g, innerpoints_g, p=2)
D = torch.gt(radius_g, D).type(torch.int)
D = D.cpu().numpy()
# Done

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
save_obj("./output/mesh_selected_inner_points.obj", inner_points[value_pos])





# Author: Frank ZY Dou
import os
import torch
import trimesh
import numpy as np
from tqdm import tqdm
from utils import  save_obj,read_VD, read_point, winding_number
from scipy.optimize import milp, Bounds, LinearConstraint

real_name = '01Ants-12_pc'
dilation = 0.025
inner_points = "random"
max_time_SCP = 1000 # in second

point_set = trimesh.load('./input/%s.obj'%real_name)
point_set = np.array(point_set.vertices)
inner_points = trimesh.load("./input/%s_random.obj"%real_name)
inner_points = np.array(inner_points.vertices)
inner_points = inner_points[np.random.choice(np.arange(len(inner_points)), 30000)] # downsample inner points to 50000.
print("The number of sampled inner candidates: ", len(inner_points))
print("The number of surface samples: ", len(point_set))
inner_points_g = torch.tensor(inner_points).cuda().double()
point_set_g = torch.tensor(point_set).cuda().double()
dist = torch.cdist(inner_points_g, point_set_g, p=2)
radius = dist.topk(1, largest=False).values
radius = radius + dilation


save_obj("./output/pc_samples.obj", point_set) # to be covered surface samples.
save_obj("./output/pc_inner_points.obj", inner_points) # candidate inner points.

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
save_obj("./output/pc_selected_inner_points.obj", inner_points[value_pos])





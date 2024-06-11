# Author: Zimeng Wang* Zhiyang Dou*

import os
import torch
import trimesh
import numpy as np
from tqdm import tqdm
from utils import save_obj, save_txt, read_VD, winding_number

def heuristic_alg(D, candidate, radius_list, reg_radius=1, reg=1, max_iter=1000, penalty='stand'):
    m, n = D.shape
    S = np.arange(m)
    A = []
    grade = []
    for i in tqdm(range(max_iter)):
        score = np.sum(D[S], axis=0).astype(float)  # summarize each col of subarray D[S]
        score = (score - np.mean(score)) / np.std(score, ddof=1)
        if len(A) > 0:
            loss = compute_min_distances(candidate, candidate[A])
            loss = (loss - np.mean(loss)) / np.std(loss, ddof=1)
            score += reg * loss
        if penalty == 'stand':
            loss_radius = 1 / radius_list
            loss_radius = (loss_radius - np.mean(loss_radius)) / np.std(loss_radius, ddof=1)
        else:
            radius_max = np.max(radius_list)
            loss_radius = 0.1 * radius_max / radius_list
        score -= reg_radius * loss_radius
        i_k = np.argmax(score)
        A.append(i_k)
        grade.append(score[i_k])
        S = S[D[S, i_k] == 0]
        if len(S) == 0:
            break
    coverage_rate = len(S) / m
    A = np.array(A)
    return A, grade, coverage_rate


def compute_min_distances(X, selected_pts):
    distances = np.linalg.norm(X[:, np.newaxis] - selected_pts, axis=2)
    min_distances = np.min(distances, axis=1)
    return min_distances


real_name = '01Ants-12_mesh'
surface_sample_num = 1500
dilation = 0.02
# inner_points = "voronoi"
inner_points = "random"

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
radius_ori = radius.cpu().numpy()
radius_list = np.reshape(radius_ori, -1)
radius = radius + dilation

save_obj("./output/pc_samples.obj", point_set) # to be covered surface samples.
save_obj("./output/pc_inner_points.obj", inner_points) # candidate inner points.

# Coverage Matrix -> GPU.
point_set_g = torch.tensor(point_set).cuda().double()
innerpoints_g = torch.tensor(inner_points).cuda().double()
radius_g = torch.tensor(radius).cuda().double()
radius_g = radius_g[:, 0]
radius_g = radius_g.unsqueeze(0).repeat(len(point_set), 1)
D = torch.cdist(point_set_g, innerpoints_g, p=2)
D = torch.gt(radius_g, D).type(torch.int)
D = D.cpu().numpy()
candidates = innerpoints_g.numpy()
# Done

# solve by heuristic algorithm
value_pos, grade, coverage_rate = heuristic_alg(D, candidates, radius_list, reg_radius=1, reg=1, max_iter=50, penalty='')
print("Coverage rate: ", 100*(1-coverage_rate), "%")
print("The number of selected inner points: ", len(value_pos))
save_obj("./output/pc_selected_inner_points.obj", inner_points[value_pos])
save_txt("./output/pc_selected_inner_points.txt", np.concatenate((inner_points[value_pos], radius_ori[value_pos]), axis=1))
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
集成QMAT和CoverageAxis的完整脚本
Author: Based on Coverage_Axis_plusplus_mesh.py
Updated: 根据QMAT使用指南完善
"""

import os
import sys
import argparse
import shutil
import subprocess
import glob
from pathlib import Path

try:
    import torch
    import trimesh
    import numpy as np
    from tqdm import tqdm
    from utils import save_obj, save_txt, read_VD, winding_number
    DEPENDENCIES_AVAILABLE = True
except ImportError as e:
    print(f"警告：缺少依赖库: {e}")
    print("请安装: pip install torch trimesh numpy tqdm")
    DEPENDENCIES_AVAILABLE = False


def heuristic_alg(D, candidate, radius_list, reg_radius=1, reg=1, max_iter=1000, penalty='stand'):
    """启发式算法求解覆盖问题"""
    m, n = D.shape
    S = np.arange(m)
    A = []
    grade = []
    pbar = tqdm(range(max_iter))
    for i in pbar:
        score = np.sum(D[S], axis=0).astype(float)
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
        pbar.set_description(f'Coverage rate: {1 - len(S) / m:.4f}')
    coverage_rate = len(S) / m
    A = np.array(A)
    return A, grade, coverage_rate


def compute_min_distances(X, selected_pts):
    """计算最小距离"""
    distances = np.linalg.norm(X[:, np.newaxis] - selected_pts, axis=2)
    min_distances = np.min(distances, axis=1)
    return min_distances


def extract_vertices_from_ma(input_file, output_file):
    """从.ma文件中提取顶点信息，保存为VD格式"""
    vertices = []
    
    with open(input_file, 'r') as f:
        lines = f.readlines()
    
    # 跳过第一行，处理顶点行
    for line in lines[1:]:
        line = line.strip()
        if line.startswith('v '):
            parts = line.split()
            if len(parts) >= 4:
                x, y, z = float(parts[1]), float(parts[2]), float(parts[3])
                if len(parts) >= 5:
                    r = float(parts[4])
                    vertices.append([x, y, z, r])
                else:
                    # 如果没有半径信息，设置默认值
                    vertices.append([x, y, z, 0.1])
    
    # 保存为VD格式
    with open(output_file, 'w') as f:
        for vertex in vertices:
            f.write(f"v {vertex[0]} {vertex[1]} {vertex[2]} {vertex[3]}\n")
    
    print(f"从 {input_file} 提取了 {len(vertices)} 个顶点，保存到 {output_file}")
    return len(vertices)


def save_selected_points_for_qmat(points_with_radius, output_file):
    """将选择的点保存为QMAT需要的格式"""
    with open(output_file, 'w') as f:
        for point in points_with_radius:
            f.write(f"v {point[0]} {point[1]} {point[2]} {point[3]}\n")
    print(f"保存了 {len(points_with_radius)} 个选择的点到 {output_file}")


def run_qmat_step1(qmat_path, input_mesh_path, input_ma_path, target_vertices=500, output_dir="./qmat_output/"):
    """运行QMAT第一步：常规简化"""
    print(f"步骤1: 使用QMAT进行常规简化到 {target_vertices} 个球...")
    
    # 确保输出目录存在且以/结尾
    os.makedirs(output_dir, exist_ok=True)
    if not output_dir.endswith('/'):
        output_dir += '/'
    
    cmd = [qmat_path, "1", input_mesh_path, input_ma_path, str(target_vertices), output_dir]
    print(f"执行命令: {' '.join(cmd)}")
    
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
        print("QMAT步骤1执行成功")
        if result.stdout:
            print("输出:", result.stdout)
        
        # 查找生成的MA文件
        ma_files = glob.glob(os.path.join(output_dir, "export_half___v_*___e_*___f_*.ma"))
        if ma_files:
            simplified_ma_file = ma_files[0]
            print(f"找到简化的MA文件: {simplified_ma_file}")
            return True, simplified_ma_file
        else:
            print("警告：未找到简化的MA文件")
            return True, None
            
    except subprocess.CalledProcessError as e:
        print(f"QMAT步骤1执行失败: {e}")
        if e.stderr:
            print("错误信息:", e.stderr)
        return False, None


def run_coverage_axis(input_mesh_path, vd_file_path, surface_sample_num=3000, dilation=0.05):
    """运行CoverageAxis算法"""
    print("步骤2: 运行CoverageAxis算法...")
    
    if not DEPENDENCIES_AVAILABLE:
        print("错误：缺少必要的依赖库，无法运行CoverageAxis")
        return False
    
    # 加载mesh
    mesh = trimesh.load(input_mesh_path)
    point_set = trimesh.sample.sample_surface(mesh, surface_sample_num)
    
    mesh_faces = np.array(mesh.faces)
    mesh_vertices = np.array(mesh.vertices)
    point_set = np.array(point_set[0])
    
    print(f"Mesh信息: 面数={mesh_faces.shape[0]}, 顶点数={mesh_vertices.shape[0]}, 采样点数={point_set.shape[0]}")
    
    # 读取VD文件
    try:
        inner_points, radius = read_VD(vd_file_path)
        inner_points = np.array(inner_points)
        radius_ori = np.array(radius)
        radius = radius_ori + dilation
        radius_list = np.reshape(radius_ori, -1)
        print(f"从VD文件读取了 {len(inner_points)} 个内部点")
    except Exception as e:
        print(f"读取VD文件失败: {e}")
        return False
    
    # 确保输出目录存在
    os.makedirs("./output", exist_ok=True)
    
    # 保存中间结果
    save_obj("./output/mesh.obj", mesh_vertices, mesh_faces)
    save_obj(f"./output/mesh_samples_{surface_sample_num}.obj", point_set)
    save_obj("./output/mesh_inner_points.obj", inner_points)
    
    # 计算覆盖矩阵
    print("计算覆盖矩阵...")
    point_set_g = torch.tensor(point_set).cuda().double()
    innerpoints_g = torch.tensor(inner_points).cuda().double()
    radius_g = torch.tensor(radius).cuda().double()
    radius_g = radius_g[:, 0]
    radius_g = radius_g.unsqueeze(0).repeat(len(point_set), 1)
    D = torch.cdist(point_set_g, innerpoints_g, p=2)
    D = torch.gt(radius_g, D).type(torch.int)
    D = D.cpu().numpy()
    candidates = innerpoints_g.cpu().numpy()
    
    # 使用启发式算法求解
    print("使用启发式算法求解覆盖问题...")
    value_pos, grade, coverage_rate = heuristic_alg(D, candidates, radius_list, 
                                                   reg_radius=1, reg=1, max_iter=100, penalty='')
    
    print(f"覆盖率: {100*(1-coverage_rate):.2f}%")
    print(f"选择的内部点数量: {len(value_pos)}")
    
    # 保存结果
    selected_points = inner_points[value_pos]
    selected_radius = radius_ori[value_pos]
    
    save_obj("./output/mesh_selected_inner_points.obj", selected_points)
    save_txt("./output/mesh_selected_inner_points.txt", 
             np.concatenate((selected_points, selected_radius), axis=1))
    
    # 为QMAT保存选择的点（格式：v x y z r）
    points_with_radius = np.concatenate((selected_points, selected_radius), axis=1)
    save_selected_points_for_qmat(points_with_radius, "./output/selected_points_for_qmat.txt")
    
    return True


def run_qmat_step2(qmat_path, input_mesh_path, input_ma_path, target_vertices, 
                   selected_points_file, output_dir="./final_output/"):
    """运行QMAT第二步：使用选择的极点进行简化"""
    print("步骤3: 使用QMAT进行带选择极点的简化...")
    
    # 确保输出目录存在且以/结尾
    os.makedirs(output_dir, exist_ok=True)
    if not output_dir.endswith('/'):
        output_dir += '/'
    
    cmd = [qmat_path, "2", input_mesh_path, input_ma_path, str(target_vertices), 
           output_dir, selected_points_file]
    print(f"执行命令: {' '.join(cmd)}")
    
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
        print("QMAT步骤2执行成功")
        if result.stdout:
            print("输出:", result.stdout)
        
        # 查找生成的文件
        obj_files = glob.glob(os.path.join(output_dir, "sim_MA___v_*___e_*___f_*.obj"))
        ma_files = glob.glob(os.path.join(output_dir, "export_half___v_*___e_*___f_*.ma"))
        poles_files = glob.glob(os.path.join(output_dir, "test_all_poles.obj"))
        
        print("生成的文件:")
        if obj_files:
            print(f"- 简化的MA (OBJ): {obj_files[0]}")
        if ma_files:
            print(f"- 简化的MA (MA): {ma_files[0]}")
        if poles_files:
            print(f"- 所有极点可视化: {poles_files[0]}")
        
        return True, obj_files[0] if obj_files else None, ma_files[0] if ma_files else None
        
    except subprocess.CalledProcessError as e:
        print(f"QMAT步骤2执行失败: {e}")
        if e.stderr:
            print("错误信息:", e.stderr)
        return False, None, None


def main():
    """主函数"""
    parser = argparse.ArgumentParser(description='集成QMAT和CoverageAxis的完整流程')
    parser.add_argument('--mesh', required=True, help='输入mesh文件路径 (.off)')
    parser.add_argument('--ma', required=True, help='输入MA文件路径 (.ma)')
    parser.add_argument('--qmat', required=True, help='QMAT可执行文件路径')
    parser.add_argument('--vertices', type=int, default=500, help='目标球数量 (默认: 500)')
    parser.add_argument('--samples', type=int, default=3000, help='表面采样点数量 (默认: 3000)')
    parser.add_argument('--dilation', type=float, default=0.05, help='膨胀参数 (默认: 0.05)')
    parser.add_argument('--temp-dir', default='./qmat_temp/', help='QMAT临时输出目录 (默认: ./qmat_temp/)')
    parser.add_argument('--output-dir', default='./final_output/', help='最终输出目录 (默认: ./final_output/)')
    parser.add_argument('--skip-step1', action='store_true', help='跳过QMAT步骤1，直接使用原始MA文件')
    
    args = parser.parse_args()
    
    # 检查输入文件
    if not os.path.exists(args.mesh):
        print(f"错误：mesh文件不存在: {args.mesh}")
        return False
    
    if not os.path.exists(args.ma):
        print(f"错误：MA文件不存在: {args.ma}")
        return False
    
    if not os.path.exists(args.qmat):
        print(f"错误：QMAT可执行文件不存在: {args.qmat}")
        return False
    
    # 创建必要的目录
    os.makedirs("./input", exist_ok=True)
    os.makedirs("./output", exist_ok=True)
    
    # 获取文件名（不含扩展名）
    mesh_name = Path(args.mesh).stem
    
    # 定义文件路径
    vd_file = f"./input/{mesh_name}_VD.txt"
    selected_points_file = "./output/selected_points_for_qmat.txt"
    
    print("="*60)
    print("集成QMAT和CoverageAxis流程开始")
    print("="*60)
    print(f"输入mesh: {args.mesh}")
    print(f"输入MA: {args.ma}")
    print(f"QMAT路径: {args.qmat}")
    print(f"目标球数量: {args.vertices}")
    print(f"临时目录: {args.temp_dir}")
    print(f"最终输出目录: {args.output_dir}")
    print("="*60)
    
    try:
        simplified_ma_file = None
        
        if not args.skip_step1:
            # 步骤1: 使用QMAT进行常规简化
            success, simplified_ma_file = run_qmat_step1(args.qmat, args.mesh, args.ma, 
                                                        args.vertices, args.temp_dir)
            if not success:
                print("步骤1失败，流程终止")
                return False
        
        # 提取VD文件
        if simplified_ma_file and os.path.exists(simplified_ma_file):
            print(f"从简化的MA文件提取VD: {simplified_ma_file}")
            extract_vertices_from_ma(simplified_ma_file, vd_file)
        else:
            print(f"从原始MA文件提取VD: {args.ma}")
            extract_vertices_from_ma(args.ma, vd_file)
        
        # 步骤2: 运行CoverageAxis
        if not run_coverage_axis(args.mesh, vd_file, args.samples, args.dilation):
            print("步骤2失败，流程终止")
            return False
        
        # 检查选择的点文件是否存在
        if not os.path.exists(selected_points_file):
            print(f"错误：未找到选择的点文件: {selected_points_file}")
            return False
        
        # 步骤3: 使用QMAT进行带选择极点的简化
        success, final_obj, final_ma = run_qmat_step2(args.qmat, args.mesh, args.ma, 
                                                     args.vertices, selected_points_file, 
                                                     args.output_dir)
        if not success:
            print("步骤3失败，流程终止")
            return False
        
        print("="*60)
        print("流程完成！")
        print("生成的文件:")
        if final_obj:
            print(f"- 最终简化MA (OBJ): {final_obj}")
        if final_ma:
            print(f"- 最终简化MA (MA): {final_ma}")
        print(f"- 中间结果目录: ./output/")
        print(f"- 最终结果目录: {args.output_dir}")
        print("="*60)
        
        return True
        
    except Exception as e:
        print(f"流程执行过程中出现错误: {e}")
        import traceback
        traceback.print_exc()
        return False


if __name__ == "__main__":
    success = main()
    sys.exit(0 if success else 1) 
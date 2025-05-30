#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Integrated QMAT and CoverageAxis complete script
Author: Based on Coverage_Axis_plusplus_mesh.py
Updated: Improved according to QMAT usage guide, supports independent output directories for multiple runs
"""

import os
import sys
import argparse
import shutil
import subprocess
import glob
from pathlib import Path
from datetime import datetime

try:
    import torch
    import trimesh
    import numpy as np
    from tqdm import tqdm
    from utils import save_obj, save_txt, read_VD, winding_number
    DEPENDENCIES_AVAILABLE = True
except ImportError as e:
    print(f"Warning: Missing dependency library: {e}")
    print("Please install: pip install torch trimesh numpy tqdm")
    DEPENDENCIES_AVAILABLE = False


def heuristic_alg(D, candidate, radius_list, reg_radius=1, reg=1, max_iter=1000, penalty='stand'):
    """Heuristic algorithm for solving coverage problem"""
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
    """Compute minimum distances"""
    distances = np.linalg.norm(X[:, np.newaxis] - selected_pts, axis=2)
    min_distances = np.min(distances, axis=1)
    return min_distances


def create_run_directory(mesh_path, base_output_dir="./runs"):
    """Create independent output directory for each run"""
    # Get mesh filename (without extension)
    mesh_name = Path(mesh_path).stem
    
    # Create timestamp
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    
    # Create run directory name
    run_dir_name = f"{mesh_name}_{timestamp}"
    run_dir = os.path.join(base_output_dir, run_dir_name)
    
    # Create directory structure
    os.makedirs(run_dir, exist_ok=True)
    os.makedirs(os.path.join(run_dir, "input"), exist_ok=True)
    os.makedirs(os.path.join(run_dir, "coverage_axis_output"), exist_ok=True)
    os.makedirs(os.path.join(run_dir, "qmat_temp"), exist_ok=True)
    os.makedirs(os.path.join(run_dir, "final_output"), exist_ok=True)
    
    print(f"Created run directory: {run_dir}")
    return run_dir


def extract_vertices_from_ma(input_file, output_file):
    """Extract vertex information from .ma file and save in VD format"""
    vertices = []
    
    with open(input_file, 'r') as f:
        lines = f.readlines()
    
    # Skip first line, process vertex lines
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
                    # If no radius information, set default value
                    vertices.append([x, y, z, 0.1])
    
    # Save in VD format
    with open(output_file, 'w') as f:
        for vertex in vertices:
            f.write(f"v {vertex[0]} {vertex[1]} {vertex[2]} {vertex[3]}\n")
    
    print(f"Extracted {len(vertices)} vertices from {input_file}, saved to {output_file}")
    return len(vertices)


def save_selected_points_for_qmat(points_with_radius, output_file):
    """Save selected points in the format required by QMAT"""
    with open(output_file, 'w') as f:
        for point in points_with_radius:
            f.write(f"v {point[0]} {point[1]} {point[2]} {point[3]}\n")
    print(f"Saved {len(points_with_radius)} selected points to {output_file}")


def run_qmat_step1(qmat_path, input_mesh_path, input_ma_path, target_vertices=500, output_dir="./qmat_output/"):
    """Run QMAT step 1: Regular simplification"""
    print(f"Step 1: Using QMAT for regular simplification to {target_vertices} spheres...")
    
    # Ensure output directory exists and ends with /
    os.makedirs(output_dir, exist_ok=True)
    if not output_dir.endswith('/'):
        output_dir += '/'
    
    cmd = [qmat_path, "1", input_mesh_path, input_ma_path, str(target_vertices), output_dir]
    print(f"Executing command: {' '.join(cmd)}")
    
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
        print("QMAT step 1 executed successfully")
        if result.stdout:
            print("Output:", result.stdout)
        
        # Find generated MA files
        ma_files = glob.glob(os.path.join(output_dir, "export_half___v_*___e_*___f_*.ma"))
        if ma_files:
            simplified_ma_file = ma_files[0]
            print(f"Found simplified MA file: {simplified_ma_file}")
            return True, simplified_ma_file
        else:
            print("Warning: Simplified MA file not found")
            return True, None
            
    except subprocess.CalledProcessError as e:
        print(f"QMAT step 1 execution failed: {e}")
        if e.stderr:
            print("Error message:", e.stderr)
        return False, None


def run_coverage_axis(input_mesh_path, vd_file_path, output_dir, surface_sample_num=3000, dilation=0.05):
    """Run CoverageAxis algorithm"""
    print("Step 2: Running CoverageAxis algorithm...")
    
    if not DEPENDENCIES_AVAILABLE:
        print("Error: Missing necessary dependency libraries, cannot run CoverageAxis")
        return False
    
    # Load mesh
    mesh = trimesh.load(input_mesh_path)
    point_set = trimesh.sample.sample_surface(mesh, surface_sample_num)
    
    mesh_faces = np.array(mesh.faces)
    mesh_vertices = np.array(mesh.vertices)
    point_set = np.array(point_set[0])
    
    print(f"Mesh info: faces={mesh_faces.shape[0]}, vertices={mesh_vertices.shape[0]}, sampling points={point_set.shape[0]}")
    
    # Read VD file
    try:
        inner_points, radius = read_VD(vd_file_path)
        inner_points = np.array(inner_points)
        radius_ori = np.array(radius)
        radius = radius_ori + dilation
        radius_list = np.reshape(radius_ori, -1)
        print(f"Read {len(inner_points)} interior points from VD file")
    except Exception as e:
        print(f"Failed to read VD file: {e}")
        return False
    
    # Ensure output directory exists
    os.makedirs(output_dir, exist_ok=True)
    
    # Save intermediate results
    save_obj(os.path.join(output_dir, "mesh.obj"), mesh_vertices, mesh_faces)
    save_obj(os.path.join(output_dir, f"mesh_samples_{surface_sample_num}.obj"), point_set)
    save_obj(os.path.join(output_dir, "mesh_inner_points.obj"), inner_points)
    
    # Calculate coverage matrix
    print("Calculating coverage matrix...")
    point_set_g = torch.tensor(point_set).cuda().double()
    innerpoints_g = torch.tensor(inner_points).cuda().double()
    radius_g = torch.tensor(radius).cuda().double()
    radius_g = radius_g[:, 0]
    radius_g = radius_g.unsqueeze(0).repeat(len(point_set), 1)
    D = torch.cdist(point_set_g, innerpoints_g, p=2)
    D = torch.gt(radius_g, D).type(torch.int)
    D = D.cpu().numpy()
    candidates = innerpoints_g.cpu().numpy()
    
    # Solve using heuristic algorithm
    print("Solving coverage problem using heuristic algorithm...")
    value_pos, grade, coverage_rate = heuristic_alg(D, candidates, radius_list, 
                                                   reg_radius=1, reg=1, max_iter=100, penalty='')
    
    print(f"Coverage rate: {100*(1-coverage_rate):.2f}%")
    print(f"Number of selected interior points: {len(value_pos)}")
    
    # Save results
    selected_points = inner_points[value_pos]
    selected_radius = radius_ori[value_pos]
    
    save_obj(os.path.join(output_dir, "mesh_selected_inner_points.obj"), selected_points)
    save_txt(os.path.join(output_dir, "mesh_selected_inner_points.txt"), 
             np.concatenate((selected_points, selected_radius), axis=1))
    
    # Save selected points for QMAT (format: v x y z r)
    points_with_radius = np.concatenate((selected_points, selected_radius), axis=1)
    selected_points_file = os.path.join(output_dir, "selected_points_for_qmat.txt")
    save_selected_points_for_qmat(points_with_radius, selected_points_file)
    
    return True, selected_points_file


def run_qmat_step2(qmat_path, input_mesh_path, input_ma_path, target_vertices, 
                   selected_points_file, output_dir="./final_output/"):
    """Run QMAT step 2: Simplification using selected poles"""
    print("Step 3: Using QMAT for simplification with selected poles...")
    
    # Ensure output directory exists and ends with /
    os.makedirs(output_dir, exist_ok=True)
    if not output_dir.endswith('/'):
        output_dir += '/'
    
    cmd = [qmat_path, "2", input_mesh_path, input_ma_path, str(target_vertices), 
           output_dir, selected_points_file]
    print(f"Executing command: {' '.join(cmd)}")
    
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
        print("QMAT step 2 executed successfully")
        if result.stdout:
            print("Output:", result.stdout)
        
        # Find generated files
        obj_files = glob.glob(os.path.join(output_dir, "sim_MA___v_*___e_*___f_*.obj"))
        ma_files = glob.glob(os.path.join(output_dir, "export_half___v_*___e_*___f_*.ma"))
        poles_files = glob.glob(os.path.join(output_dir, "test_all_poles.obj"))
        
        print("Generated files:")
        if obj_files:
            print(f"- Simplified MA (OBJ): {obj_files[0]}")
        if ma_files:
            print(f"- Simplified MA (MA): {ma_files[0]}")
        if poles_files:
            print(f"- All poles visualization: {poles_files[0]}")
        
        return True, obj_files[0] if obj_files else None, ma_files[0] if ma_files else None
        
    except subprocess.CalledProcessError as e:
        print(f"QMAT step 2 execution failed: {e}")
        if e.stderr:
            print("Error message:", e.stderr)
        return False, None, None


def save_run_info(run_dir, args, results):
    """Save run information to file"""
    info_file = os.path.join(run_dir, "run_info.txt")
    with open(info_file, 'w', encoding='utf-8') as f:
        f.write("="*60 + "\n")
        f.write("Run Information\n")
        f.write("="*60 + "\n")
        f.write(f"Run time: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
        f.write(f"Input mesh: {args.mesh}\n")
        f.write(f"Input MA: {args.ma}\n")
        f.write(f"QMAT path: {args.qmat}\n")
        f.write(f"Target number of spheres: {args.vertices}\n")
        f.write(f"Surface sampling points: {args.samples}\n")
        f.write(f"Dilation parameter: {args.dilation}\n")
        f.write(f"Skip step 1: {args.skip_step1}\n")
        f.write("\n")
        f.write("Generated files:\n")
        for key, value in results.items():
            if value:
                f.write(f"- {key}: {value}\n")
        f.write("="*60 + "\n")
    
    print(f"Run information saved to: {info_file}")


def main():
    """Main function"""
    parser = argparse.ArgumentParser(description='Complete pipeline integrating QMAT and CoverageAxis')
    parser.add_argument('--mesh', required=True, help='Input mesh file path (.off)')
    parser.add_argument('--ma', required=True, help='Input MA file path (.ma)')
    parser.add_argument('--qmat', required=True, help='QMAT executable file path')
    parser.add_argument('--vertices', type=int, default=500, help='Target number of spheres (default: 500)')
    parser.add_argument('--samples', type=int, default=3000, help='Surface sampling points (default: 3000)')
    parser.add_argument('--dilation', type=float, default=0.05, help='Dilation parameter (default: 0.05)')
    parser.add_argument('--runs-dir', default='./runs', help='Run directory root (default: ./runs)')
    parser.add_argument('--skip-step1', action='store_true', help='Skip QMAT step 1, directly use original MA file')
    
    args = parser.parse_args()
    
    # Check input files
    if not os.path.exists(args.mesh):
        print(f"Error: mesh file does not exist: {args.mesh}")
        return False
    
    if not os.path.exists(args.ma):
        print(f"Error: MA file does not exist: {args.ma}")
        return False
    
    if not os.path.exists(args.qmat):
        print(f"Error: QMAT executable does not exist: {args.qmat}")
        return False
    
    # Create run directory
    run_dir = create_run_directory(args.mesh, args.runs_dir)
    
    # Define subdirectories
    input_dir = os.path.join(run_dir, "input")
    coverage_output_dir = os.path.join(run_dir, "coverage_axis_output")
    qmat_temp_dir = os.path.join(run_dir, "qmat_temp")
    final_output_dir = os.path.join(run_dir, "final_output")
    
    # Get filename (without extension)
    mesh_name = Path(args.mesh).stem
    
    # Define file paths
    vd_file = os.path.join(input_dir, f"{mesh_name}_VD.txt")
    
    print("="*60)
    print("Integrated QMAT and CoverageAxis pipeline started")
    print("="*60)
    print(f"Input mesh: {args.mesh}")
    print(f"Input MA: {args.ma}")
    print(f"QMAT path: {args.qmat}")
    print(f"Target number of spheres: {args.vertices}")
    print(f"Run directory: {run_dir}")
    print("="*60)
    
    # For saving result information
    results = {}
    
    try:
        simplified_ma_file = None
        
        if not args.skip_step1:
            # Step 1: Use QMAT for regular simplification
            success, simplified_ma_file = run_qmat_step1(args.qmat, args.mesh, args.ma, 
                                                        args.vertices, qmat_temp_dir + "/")
            if not success:
                print("Step 1 failed, pipeline terminated")
                return False
            results["QMAT step 1 simplified MA"] = simplified_ma_file
        
        # Extract VD file
        if simplified_ma_file and os.path.exists(simplified_ma_file):
            print(f"Extracting VD from simplified MA file: {simplified_ma_file}")
            extract_vertices_from_ma(simplified_ma_file, vd_file)
        else:
            print(f"Extracting VD from original MA file: {args.ma}")
            extract_vertices_from_ma(args.ma, vd_file)
        results["VD file"] = vd_file
        
        # Step 2: Run CoverageAxis
        coverage_result = run_coverage_axis(args.mesh, vd_file, coverage_output_dir, 
                                          args.samples, args.dilation)
        if not coverage_result:
            print("Step 2 failed, pipeline terminated")
            return False
        
        success, selected_points_file = coverage_result
        results["Selected points file"] = selected_points_file
        
        # Check if selected points file exists
        if not os.path.exists(selected_points_file):
            print(f"Error: Selected points file not found: {selected_points_file}")
            return False
        
        # Step 3: Use QMAT for simplification with selected poles
        success, final_obj, final_ma = run_qmat_step2(args.qmat, args.mesh, args.ma, 
                                                     args.vertices, selected_points_file, 
                                                     final_output_dir + "/")
        if not success:
            print("Step 3 failed, pipeline terminated")
            return False
        
        results["Final simplified MA (OBJ)"] = final_obj
        results["Final simplified MA (MA)"] = final_ma
        
        # Save run information
        save_run_info(run_dir, args, results)
        
        print("="*60)
        print("Pipeline completed!")
        print(f"Run directory: {run_dir}")
        print("Directory structure:")
        print(f"├── input/                    # Input and VD files")
        print(f"├── coverage_axis_output/     # CoverageAxis intermediate results")
        print(f"├── qmat_temp/               # QMAT step 1 temporary files")
        print(f"├── final_output/            # Final output files")
        print(f"└── run_info.txt             # Run information record")
        print("="*60)
        
        return True
        
    except Exception as e:
        print(f"Error occurred during pipeline execution: {e}")
        import traceback
        traceback.print_exc()
        return False


if __name__ == "__main__":
    success = main()
    sys.exit(0 if success else 1) 
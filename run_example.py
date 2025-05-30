#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Usage example script
Updated: Updated according to QMAT usage guide
"""

import subprocess
import sys
import os

def run_integrated_pipeline():
    """Example of running integrated pipeline"""
    
    # Configuration file paths
    mesh_file = "./input/offs/exp-lbs-beagle.off"
    ma_file = "./input/offs/exp-lbs-beagle.ma"            # Input MA file
    qmat_executable = "./skel_connection/QMAT/build/QMAT"  # QMAT executable path
    
    # Check if files exist
    if not os.path.exists(mesh_file):
        print(f"Error: mesh file does not exist: {mesh_file}")
        print("Please ensure input file is in .off format")
        return False
    
    if not os.path.exists(ma_file):
        print(f"Error: MA file does not exist: {ma_file}")
        return False
    
    if not os.path.exists(qmat_executable):
        print(f"Error: QMAT executable does not exist: {qmat_executable}")
        return False
    
    # Build command
    cmd = [
        "python", "integrated_qmat_coverage_axis.py",
        "--mesh", mesh_file,
        "--ma", ma_file,
        "--qmat", qmat_executable,
        "--vertices", "500",           # Target number of spheres
        "--samples", "3000",           # Surface sampling points
        "--dilation", "0.05",          # Dilation parameter
    ]
    
    print("Running integrated pipeline...")
    print(f"Command: {' '.join(cmd)}")
    print("\nPipeline description:")
    print("1. QMAT regular simplification -> Generate simplified MA file")
    print("2. CoverageAxis algorithm -> Select optimal interior points")
    print("3. QMAT simplification with selected poles -> Generate final result")
    print()
    
    try:
        result = subprocess.run(cmd, check=True)
        print("\nPipeline execution successful!")
        print("\nGenerated files:")
        print("- Intermediate results: ./output/")
        print("- QMAT temporary files: ./qmat_temp/")
        print("- Final results: ./final_output/")
        return True
    except subprocess.CalledProcessError as e:
        print(f"\nPipeline execution failed: {e}")
        return False


def run_skip_step1_example():
    """Example of skipping QMAT step 1 (directly use original MA file)"""
    
    mesh_file = "./input/offs/exp-lbs-beagle.off"
    ma_file = "./input/offs/exp-lbs-beagle.ma"
    qmat_executable = "./skel_connection/QMAT/build/QMAT"
    
    cmd = [
        "python", "integrated_qmat_coverage_axis.py",
        "--mesh", mesh_file,
        "--ma", ma_file,
        "--qmat", qmat_executable,
        "--vertices", "500",
    ]
    
    print("Running pipeline that skips step 1...")
    print(f"Command: {' '.join(cmd)}")
    
    try:
        result = subprocess.run(cmd, check=True)
        print("Pipeline execution successful!")
        return True
    except subprocess.CalledProcessError as e:
        print(f"Pipeline execution failed: {e}")
        return False


if __name__ == "__main__":
    print("Select running mode:")
    print("1. Complete pipeline (including QMAT step 1)")
    print("2. Pipeline skipping QMAT step 1")
    
    choice = input("Please enter your choice (1 or 2): ").strip()
    
    if choice == "1":
        success = run_integrated_pipeline()
    elif choice == "2":
        success = run_skip_step1_example()
    else:
        print("Invalid choice, running default complete pipeline")
        success = run_integrated_pipeline()
    
    sys.exit(0 if success else 1) 
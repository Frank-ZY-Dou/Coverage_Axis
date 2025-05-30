#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
使用示例脚本
Updated: 根据QMAT使用指南更新
"""

import subprocess
import sys
import os

def run_integrated_pipeline():
    """运行集成流程的示例"""
    
    # 配置文件路径
    mesh_file = "./input/offs/exp-lbs-beagle.off"
    ma_file = "./input/offs/exp-lbs-beagle.ma"            # 输入MA文件
    qmat_executable = "./skel_connection/QMAT/build/QMAT"  # QMAT可执行文件路径
    
    # 检查文件是否存在
    if not os.path.exists(mesh_file):
        print(f"错误：mesh文件不存在: {mesh_file}")
        print("请确保输入文件为.off格式")
        return False
    
    if not os.path.exists(ma_file):
        print(f"错误：MA文件不存在: {ma_file}")
        return False
    
    if not os.path.exists(qmat_executable):
        print(f"错误：QMAT可执行文件不存在: {qmat_executable}")
        return False
    
    # 构建命令
    cmd = [
        "python", "integrated_qmat_coverage_axis.py",
        "--mesh", mesh_file,
        "--ma", ma_file,
        "--qmat", qmat_executable,
        "--vertices", "500",           # 目标球数量
        "--samples", "3000",           # 表面采样点数量
        "--dilation", "0.05",          # 膨胀参数
    ]
    
    print("运行集成流程...")
    print(f"命令: {' '.join(cmd)}")
    print("\n流程说明:")
    print("1. QMAT常规简化 -> 生成简化的MA文件")
    print("2. CoverageAxis算法 -> 选择最优内部点")
    print("3. QMAT带选择极点简化 -> 生成最终结果")
    print()
    
    try:
        result = subprocess.run(cmd, check=True)
        print("\n流程执行成功！")
        print("\n生成的文件:")
        print("- 中间结果: ./output/")
        print("- QMAT临时文件: ./qmat_temp/")
        print("- 最终结果: ./final_output/")
        return True
    except subprocess.CalledProcessError as e:
        print(f"\n流程执行失败: {e}")
        return False


def run_skip_step1_example():
    """跳过QMAT步骤1的示例（直接使用原始MA文件）"""
    
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
    
    print("运行跳过步骤1的流程...")
    print(f"命令: {' '.join(cmd)}")
    
    try:
        result = subprocess.run(cmd, check=True)
        print("流程执行成功！")
        return True
    except subprocess.CalledProcessError as e:
        print(f"流程执行失败: {e}")
        return False


if __name__ == "__main__":
    print("选择运行模式:")
    print("1. 完整流程（包含QMAT步骤1）")
    print("2. 跳过QMAT步骤1的流程")
    
    choice = input("请输入选择 (1 或 2): ").strip()
    
    if choice == "1":
        success = run_integrated_pipeline()
    elif choice == "2":
        success = run_skip_step1_example()
    else:
        print("无效选择，运行默认完整流程")
        success = run_integrated_pipeline()
    
    sys.exit(0 if success else 1) 
# 如何使用骨骼连接

**Language / 语言**: [English](README_EN.md) | **中文**

`integrated_qmat_coverage_axis.py`这个脚本集成了QMAT和CoverageAxis算法，提供了一个完整的medial axis简化流程。

## 功能特点

1. **步骤1**: 使用QMAT进行常规简化（可选）
2. **步骤2**: 从MA文件提取VD信息，运行CoverageAxis算法选择最优内部点
3. **步骤3**: 使用选择的极点运行QMAT进行带约束的简化

## 依赖库

```bash
pip install torch trimesh numpy tqdm
```
确保你安装了本文件夹下的QMAT

## 文件结构

```
.
├── integrated_qmat_coverage_axis.py  # 主脚本
├── run_example.py                    # 使用示例
├── input/                            # 输入文件目录
├── output/                           # CoverageAxis中间结果目录
├── qmat_temp/                        # QMAT步骤1临时输出目录
└── final_output/                     # 最终输出目录
```

## 使用方法

### 方法1: 命令行使用

```bash
python integrated_qmat_coverage_axis.py \
    --mesh ./input/bird/bird.off \
    --ma ./input/bird/bird.ma \
    --qmat ./build/QMAT \
    --vertices 500 \
    --samples 3000 \
    --dilation 0.05 \
    --temp-dir ./qmat_temp/ \
    --output-dir ./final_output/
```

### 方法2: 跳过QMAT步骤1

如果您已经有合适的MA文件，可以跳过第一步：

```bash
python integrated_qmat_coverage_axis.py \
    --mesh ./input/bird/bird.off \
    --ma ./input/bird/bird.ma \
    --qmat ./build/QMAT \
    --vertices 500 \
    --skip-step1 \
    --output-dir ./final_output/
```

### 方法3: 使用示例脚本

1. 修改 `run_example.py` 中的文件路径
2. 运行：
```bash
python run_example.py
```

## 参数说明

- `--mesh`: 输入mesh文件路径 (必须是.off格式)
- `--ma`: 输入MA文件路径 (.ma格式)
- `--qmat`: QMAT可执行文件路径
- `--vertices`: 目标球数量 (默认: 500)
- `--samples`: 表面采样点数量 (默认: 3000)
- `--dilation`: 膨胀参数 (默认: 0.05)
- `--temp-dir`: QMAT临时输出目录 (默认: ./qmat_temp/)
- `--output-dir`: 最终输出目录 (默认: ./final_output/)
- `--skip-step1`: 跳过QMAT步骤1，直接使用原始MA文件

## 输出文件

### CoverageAxis中间文件 (./output/)
- `mesh.obj`: 原始mesh
- `mesh_samples_N.obj`: 表面采样点
- `mesh_inner_points.obj`: 所有内部候选点
- `mesh_selected_inner_points.obj`: 选择的最优内部点
- `mesh_selected_inner_points.txt`: 选择的点的坐标和半径信息
- `selected_points_for_qmat.txt`: 为QMAT格式化的选择点文件

### QMAT临时文件 (./qmat_temp/)
- `export_half___v_X___e_Y___f_Z.ma`: 步骤1生成的简化MA文件
- `sim_MA___v_X___e_Y___f_Z.obj`: 步骤1生成的简化MA (OBJ格式)

### 最终输出文件 (./final_output/)
- `sim_MA___v_X___e_Y___f_Z.obj`: 最终简化的medial axis (OBJ格式)
- `export_half___v_X___e_Y___f_Z.ma`: 最终简化的medial axis (MA格式)
- `test_all_poles.obj`: 所有极点的可视化文件



## QMAT输出文件格式

根据QMAT的输出规范：
- `sim_MA___v_X___e_Y___f_Z.obj`: 简化的medial axis (OBJ格式)
- `export_half___v_X___e_Y___f_Z.ma`: 简化的medial axis (MA格式)
- `test_all_poles.obj`: (仅模式2) 所有极点可视化文件

其中X、Y、Z分别表示顶点数、边数、面数。




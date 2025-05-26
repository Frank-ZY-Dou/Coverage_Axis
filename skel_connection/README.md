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

## 流程详解

### 步骤1: QMAT常规简化 (可选)
```bash
./QMAT 1 <surface_mesh.off> <medial_mesh.ma> <num_target_spheres> <output_path>
```
- 对原始MA进行常规简化
- 生成简化的MA文件用于后续处理

### 步骤2: CoverageAxis算法
1. 从MA文件（简化的或原始的）提取VD信息
2. 加载surface mesh并进行表面采样
3. 计算覆盖矩阵
4. 使用启发式算法选择最优内部点
5. 将选择的点保存为QMAT需要的格式

### 步骤3: QMAT带选择极点的简化
```bash
./QMAT 2 <surface_mesh.off> <medial_mesh.ma> <num_target_spheres> <output_path> <selected_points.txt>
```
- 使用CoverageAxis选择的点作为约束
- 生成保持重要特征的简化MA

## QMAT输出文件格式

根据QMAT的输出规范：
- `sim_MA___v_X___e_Y___f_Z.obj`: 简化的medial axis (OBJ格式)
- `export_half___v_X___e_Y___f_Z.ma`: 简化的medial axis (MA格式)
- `test_all_poles.obj`: (仅模式2) 所有极点可视化文件

其中X、Y、Z分别表示顶点数、边数、面数。

## 注意事项

1. **输入格式要求**: 
   - Surface mesh必须是.off格式
   - Medial axis必须是.ma格式
2. **GPU要求**: CoverageAxis部分需要CUDA支持的GPU
3. **内存要求**: 大型mesh可能需要大量内存
4. **QMAT路径**: 确保QMAT可执行文件路径正确且有执行权限
5. **输出目录**: 所有输出目录路径必须以'/'结尾


## 自定义配置

可以通过修改脚本中的参数来调整算法行为：

- `surface_sample_num`: 表面采样点数量
- `dilation`: 膨胀参数
- `max_iter`: 启发式算法最大迭代次数
- `reg_radius`, `reg`: 正则化参数


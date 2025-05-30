# How to use skeleton connection

**Language / 语言**: **English** | [中文](README.md)

This script integrates QMAT and CoverageAxis algorithms, providing a complete medial axis simplification pipeline.

## Features

1. **Step 1**: Use QMAT for regular simplification (optional)
2. **Step 2**: Extract VD information from MA file, run CoverageAxis algorithm to select optimal interior points
3. **Step 3**: Use selected poles to run QMAT for constrained simplification

## Dependencies

```bash
pip install torch trimesh numpy tqdm
```
Ensure that you have installed the QMAT in this folder

## File Structure

```
.
├── integrated_qmat_coverage_axis.py  # Main script
├── run_example.py                    # Usage example
├── utils.py                          # Utility functions (copy from original project)
├── input/                            # Input files directory
├── output/                           # CoverageAxis intermediate results directory
├── qmat_temp/                        # QMAT step 1 temporary output directory
└── final_output/                     # Final output directory
```

## Usage

### Method 1: Command Line Usage

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

### Method 2: Skip QMAT Step 1

If you already have a suitable MA file, you can skip the first step:

```bash
python integrated_qmat_coverage_axis.py \
    --mesh ./input/bird/bird.off \
    --ma ./input/bird/bird.ma \
    --qmat ./build/QMAT \
    --vertices 500 \
    --skip-step1 \
    --output-dir ./final_output/
```

### Method 3: Using Example Script

1. Modify file paths in `run_example.py`
2. Run:
```bash
python run_example.py
```

## Parameter Description

- `--mesh`: Input mesh file path (must be .off format)
- `--ma`: Input MA file path (.ma format)
- `--qmat`: QMAT executable file path
- `--vertices`: Target number of spheres (default: 500)
- `--samples`: Number of surface sampling points (default: 3000)
- `--dilation`: Dilation parameter (default: 0.05)
- `--temp-dir`: QMAT temporary output directory (default: ./qmat_temp/)
- `--output-dir`: Final output directory (default: ./final_output/)
- `--skip-step1`: Skip QMAT step 1, use original MA file directly

## Output Files

### CoverageAxis Intermediate Files (./output/)
- `mesh.obj`: Original mesh
- `mesh_samples_N.obj`: Surface sampling points
- `mesh_inner_points.obj`: All interior candidate points
- `mesh_selected_inner_points.obj`: Selected optimal interior points
- `mesh_selected_inner_points.txt`: Coordinates and radius information of selected points
- `selected_points_for_qmat.txt`: Selected points file formatted for QMAT

### QMAT Temporary Files (./qmat_temp/)
- `export_half___v_X___e_Y___f_Z.ma`: Simplified MA file generated in step 1
- `sim_MA___v_X___e_Y___f_Z.obj`: Simplified MA generated in step 1 (OBJ format)

### Final Output Files (./final_output/)
- `sim_MA___v_X___e_Y___f_Z.obj`: Final simplified medial axis (OBJ format)
- `export_half___v_X___e_Y___f_Z.ma`: Final simplified medial axis (MA format)
- `test_all_poles.obj`: Visualization file of all poles


## QMAT Output File Format

According to QMAT output specifications:
- `sim_MA___v_X___e_Y___f_Z.obj`: Simplified medial axis (OBJ format)
- `export_half___v_X___e_Y___f_Z.ma`: Simplified medial axis (MA format)
- `test_all_poles.obj`: (Mode 2 only) All poles visualization file

Where X, Y, Z represent the number of vertices, edges, and faces respectively.


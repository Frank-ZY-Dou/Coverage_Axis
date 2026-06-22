# Coverage Axis SCP — GPU-accelerated exact solver

GPU-accelerated, **exact** solver for the point-selection step of
[Coverage Axis](https://github.com/Frank-ZY-Dou/Coverage_Axis) (Dou et al. 2022).

Coverage Axis picks inner skeletal points by solving a 0-1 **Set Cover Program
(SCP)** over a finite surface sample:

```
min  sum_i x_i
s.t. sum_{i : ball i covers surface point j} x_i >= 1   for every j
     x_i in {0,1}
```

where `x_i = 1` keeps inner ball `i`. This solver keeps the **same program and
its global optimum** — no heuristic, no approximation — and reaches it faster by

- building the coverage relation as packed bitsets on the GPU (NVIDIA Warp),
- solving the integer program exactly with CP-SAT (OR-Tools) or HiGHS,
- optional optimum-preserving presolve (unique-cover forcing, dominance, dedup).

## Why

The reference Coverage Axis builds the coverage matrix `D` and solves it with
`scipy.optimize.milp` (HiGHS) under a wall-clock cutoff. On larger instances
that solve becomes the bottleneck. Here the coverage matrix is built on the GPU
and the exact solve uses a parallel CP-SAT backend; the returned selection is the
**same proven-optimal** set, several times faster.

## Results (beagle surface; one RTX 2080 Ti, 8-core CP-SAT)

```
HARD instance: 8000 surface points x 6000 candidates
  coverage build : numpy 2.28s  ->  GPU/Warp 0.02s   (113x, identical coverage)
  exact solve    : HiGHS 11.66s ->  CP-SAT  3.00s    (3.9x, objective 42, proven optimal)
  end-to-end     : 13.93s       ->  3.02s            (4.6x, same global optimum)

EASY instance: 2500 x 2500   ->  objective 22, end-to-end ~1.8x, same optimum
```

All backends return the **identical objective**, proven optimal and feasible.
Notes on the numbers: CP-SAT runs with 8 workers while `scipy.milp` (HiGHS) is
serial, so part of the solve speedup is parallelism; the GPU-build and end-to-end
figures are steady-state (a single cold run pays a one-time ~7 s Warp/CUDA
startup). The selection itself is unchanged.

## Optimality is certified, not asserted

`verify.py` confirms the result independently of any solver's "optimal" flag:

- **feasibility** — recompute the coverage from raw geometry; the selection
  covers every surface point;
- **two independent solves** (HiGHS and CP-SAT) agree on the objective;
- **LP lower-bound certificate** — `ceil(LP relaxation) == objective` proves
  optimality;
- **brute force** — on small random instances the result matches the exhaustive
  minimum cover.

## Install

```bash
pip install warp-lang ortools scipy numpy        # CP-SAT backend (recommended)
# scipy alone is enough for the HiGHS backend
```

A CUDA GPU is used for the coverage build (NVIDIA Warp); a CPU fallback is
included.

## Usage

Run from this directory (uses the repo's sample mesh in `../input/`):

```bash
# 1. build an SCP instance (surface samples, inner candidates, radii) from a mesh
python gen_instance.py --off ../input/01Ants-12_mesh.off --out ants.npz \
    --n-candidates 3000 --n-surface 3000 --dilation 0.03

# 2. compare backends (solution + runtime)
python compare.py --inst ants.npz

# 3. independent optimality / feasibility check
python verify.py  --inst ants.npz --brute
```

Programmatic use:

```python
import numpy as np, scp_fast as SF
inst = np.load("ants.npz")
res = SF.coverage_axis_fast(inst["surface"], inst["centers"], inst["radii"], backend="cpsat")
selected = res["selected"]          # indices of the chosen inner points
print(res["objective"], res["proven_optimal"])
```

## Files

| file | role |
|---|---|
| `scp_fast.py`  | GPU/Warp bitset coverage build, vectorized presolve, exact-solve pipeline |
| `scp_exact.py` | coverage CSR, full presolve (dominance reductions), CP-SAT / HiGHS solvers |
| `gen_instance.py` | build an SCP instance from a `.off` surface mesh |
| `compare.py`   | head-to-head solution + runtime comparison of the backends |
| `verify.py`    | independent feasibility / LP-certificate / brute-force optimality checks |
| `bench_scp.py` | reference-vs-accelerated benchmark asserting identical optimum |

## Citation

```bibtex
@article{dou2022coverage,
  author  = {Dou, Zhiyang and Lin, Cheng and Xu, Rui and Yang, Lei and Xin, Shiqing and Komura, Taku and Wang, Wenping},
  title   = {Coverage Axis: Inner Point Selection for {3D} Shape Skeletonization},
  journal = {Computer Graphics Forum},
  volume  = {41},
  number  = {2},
  pages   = {419--432},
  year    = {2022}
}
```

## License

Released under the [MIT License](LICENSE).

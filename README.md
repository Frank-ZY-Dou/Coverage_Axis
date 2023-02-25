# Introduction
Official code for the paper [Coverage Axis: Inner Point Selection for 3D Shape Skeletonization
](https://arxiv.org/abs/2110.00965), Eurographics 2022.



Authors: [Zhiyang Dou](https://frank-zy-dou.github.io/), 
[Cheng Lin](https://clinplayer.github.io/), 
[Rui Xu](https://xrvitd.github.io/index.html), 
[Lei Yang](https://www.linkedin.cn/incareer/in/lei-yang-842052119),
[Shiqing Xin](http://irc.cs.sdu.edu.cn/~shiqing/index.html),
[Taku Komura](https://i.cs.hku.hk/~taku/), 
[Wenping Wang](https://engineering.tamu.edu/cse/profiles/Wang-Wenping.html).

![teasar](./assets/fig_teaser.jpg)
> In this paper, we present a simple yet effective formulation called Coverage Axis for 3D shape skeletonization. Inspired by the set cover problem, our key idea is to cover all the surface points using as few inside medial balls as possible. This formulation inherently induces a compact and expressive approximation of the Medial Axis Transform (MAT) of a given shape. Different from previous methods that rely on local approximation error, our method allows a global consideration of the overall shape structure, leading to an efficient high-level abstraction and superior robustness to noise. Another appealing aspect of our method is its capability to handle more generalized input such as point clouds and poor-quality meshes. Extensive comparisons and evaluations demonstrate the remarkable effectiveness of our method for generating compact and expressive skeletal representation to approximate the MAT.


This repo contains the code for **skeletal point selection** by solving SCP.

# Requirements
## System requirements
- Linux (tested on Ubuntu 20.04)
- Python 3.8
- Nvidia 3090 (GPU is used for acceleration)
## Installation

```angular2html
conda env create -f ca.yml
conda activate CA
pip install -r requirements.txt
```

# Usage


## Mesh
The input mesh `01Ants-12.off` is placed in the folder `input`. The mesh is normalized.

Specify the settings for Coverage Axis in ```Coverage_Axis_mesh.py```
```angular2html
real_name = '01Ants-12'
surface_sample_num = 2000
dilation = 0.02
inner_points = "voronoi"
```
Run
```angular2html
python Coverage_Axis_mesh.py
```
The output is placed in the folder `output`.
- `inner_points.obj` is the candidate inner points.
- `mesh.obj` is the input mesh.
- `mesh_samples_2000.obj` is the sampled surface points being covered.
- `selected_inner_points.obj` is the selected inner points.


<img src="./assets/fig_results.png" 
        alt="Picture" 
        width="800" 
        height="600" 
        style="display: block; margin: 0 auto" />

You may use randomly generated points inside the volume as inner candidate points by setting `inner_points = "random"
`. 

```angular2html
## Point Cloud

We use [Fast Winding Number](https://www.dgp.toronto.edu/projects/fast-winding-numbers/) for Inside-outside determination for point cloud inputs.
Please follow this [tutorial](https://libigl.github.io/tutorial/) for building libigl.
Note that more dependencies are needed for building libigl:
```angular2html
sudo apt-get install git
sudo apt-get install build-essential
sudo apt-get install cmake
sudo apt-get install libx11-dev
sudo apt-get install mesa-common-dev libgl1-mesa-dev libglu1-mesa-dev
sudo apt-get install libxrandr-dev
sudo apt-get install libxi-dev
sudo apt-get install libxmu-dev
sudo apt-get install libblas-dev
sudo apt-get install libxinerama-dev
sudo apt-get install libxcursor-dev
sudo apt install libeigen3-dev
sudo apt-get install libcgal-dev


git clone https://github.com/libigl/libigl.git
cd libigl/
mkdir build
cd build
cmake ../
make -j4
```

You need to run cpp codes, you also needs to install [Eigen](http://eigen.tuxfamily.org/index.php?title=Main_Page#Download). 


### Candidate Generation

We generate inside candidates based on Fast Winding Number. The candidates are generated by randomly sampling inside the volume.
Other sampling strategies like Voronoi-based sampling can also be used.

```angular2html
```







# More Information
### Solver for Coverage Axis in MATLAB
The original optimization is solved by MATLAB. In this repo, we solve SCP by Scipy in Python.
```angular2html
f =  ones(1,medial_num); 
A =  -D;
b =  -ones(boundary_num,1)*1;%here ,we fix p_i 
lb = zeros(medial_num,1);
ub = ones(medial_num,1);
iint = [1:medial_num];
tic;
[x,fval]=intlinprog(f,iint,A,b,[],[],lb,ub);
toc;
disp('min_number:');
disp(fval);
```
# References
- https://libigl.github.io/tutorial/  Many thanks to the contributors to libigl :)
- https://github.com/mayorx/hungarian-algorithm



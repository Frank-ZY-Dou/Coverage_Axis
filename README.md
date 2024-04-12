# Introduction
- 🌟Coverage Axis.
  - We provide Coverage Axis computation for both mesh and point cloud inputs. 
    - The point cloud can be **unoriented** -> check out our latest SIGGRAPH 2023 work [here](https://xrvitd.github.io/Projects/GCNO/index.html).
  - Operations are accelerated by GPU, e.g., computation of coverage matrix and winding number for a mesh.
  - We provide codes for building connectivity in [skel_connection](skel_connection): [README.md](skel_connection%2Freadme.md).
- 🌟Coverage Axis++.
  - **[New!] Coverage Axis++: Efficient Inner Point Selection for 3D Shape Skeletonization: https://arxiv.org/abs/2401.12946** 
  - Code to be release.

Some geometry tools for MAT and related topics: [Geometry_Tools](https://github.com/Frank-ZY-Dou/Geometry_Tools).

🐱 **[Coverage Axis: Inner Point Selection for 3D Shape Skeletonization
](https://arxiv.org/abs/2110.00965), Eurographics 2022**.

Authors: [Zhiyang Dou](https://frank-zy-dou.github.io/), 
[Cheng Lin](https://clinplayer.github.io/), 
[Rui Xu](https://xrvitd.github.io/index.html), 
[Lei Yang](https://www.linkedin.cn/incareer/in/lei-yang-842052119),
[Shiqing Xin](http://irc.cs.sdu.edu.cn/~shiqing/index.html),
[Taku Komura](https://i.cs.hku.hk/~taku/), 
[Wenping Wang](https://engineering.tamu.edu/cse/profiles/Wang-Wenping.html).

[[Project Page](https://frank-zy-dou.github.io/projects/CoverageAxis/index.html)][[Paper](https://arxiv.org/abs/2110.00965)][[Code](https://github.com/Frank-ZY-Dou/Coverage_Axis)]

![teasar](./assets/fig_teaser.jpg)
In this paper, we present a simple yet effective formulation called Coverage Axis for 3D shape skeletonization. Inspired by the set cover problem, our key idea is to cover all the surface points using as few inside medial balls as possible. This formulation inherently induces a compact and expressive approximation of the Medial Axis Transform (MAT) of a given shape. Different from previous methods that rely on local approximation error, our method allows a global consideration of the overall shape structure, leading to an efficient high-level abstraction and superior robustness to noise. Another appealing aspect of our method is its capability to handle more generalized input such as point clouds and poor-quality meshes. Extensive comparisons and evaluations demonstrate the remarkable effectiveness of our method for generating compact and expressive skeletal representation to approximate the MAT.

🐱 **[Coverage Axis++: Efficient Inner Point Selection for 3D Shape Skeletonization
](https://arxiv.org/abs/2401.12946), Arxiv 2023**.

Authors: Zimeng Wang*, 
[Zhiyang Dou*](https://clinplayer.github.io/),
[Rui Xu](https://xrvitd.github.io/index.html),
[Cheng Lin](https://clinplayer.github.io/), 
[Yuan Liu](https://liuyuan-pal.github.io/), 
[Xiaoxiao Long](https://www.xxlong.site/), 
[Shiqing Xin](http://irc.cs.sdu.edu.cn/~shiqing/index.html), 
[Lingjie Liu](https://lingjie0206.github.io/), 
[Taku Komura](https://www.cs.hku.hk/index.php/people/academic-staff/taku), 
[Xiaoming Yuan](https://hkumath.hku.hk/~xmyuan/),
[Wenping Wang](https://engineering.tamu.edu/cse/profiles/Wang-Wenping.html).

[[Project Page](https://frank-zy-dou.github.io/projects/CoverageAxis++/index.html)][[Paper](https://arxiv.org/pdf/2401.12946.pdf)][[Code](https://github.com/Frank-ZY-Dou/Coverage_Axis)]


![teasar](./assets/fig_teaser_cat++.png)
We introduce Coverage Axis++, a novel and efficient approach to 3D shape skeletonization. The current state-of-the-art approaches for this task often rely on the watertightness of the input or suffer from substantial computational costs, thereby limiting their practicality. To address this challenge, Coverage Axis++ proposes a heuristic algorithm to select skeletal points, offering a high-accuracy approximation of the Medial Axis Transform (MAT) while significantly mitigating computational intensity for various shape representations. We introduce a simple yet effective strategy that considers both shape coverage and uniformity to derive skeletal points. The selection procedure enforces consistency with the shape structure while favoring the dominant medial balls, which thus introduces a compact underlying shape representation in terms of MAT. As a result, Coverage Axis++ allows for skeletonization for various shape representations (e.g., water-tight meshes, triangle soups, point clouds), specification of the number of skeletal points, few hyperparameters, and highly efficient computation with improved reconstruction accuracy. Extensive experiments across a wide range of 3D shapes validate the efficiency and effectiveness of Coverage Axis++. 




# Requirements
## System requirements
- Linux Ubuntu 20.04
- Python 3.8
- Nvidia GeForce RTX 3090 (GPU is used for acceleration)
## Installation

```angular2html
conda env create -f ca.yml
conda activate CA
pip install -r requirements.txt
```

# Usage


## Mesh Input
The input mesh `01Ants-12.off` is placed in the folder `input`. The mesh is normalized.

Specify the settings for Coverage Axis in ```Coverage_Axis_mesh.py```
```angular2html
real_name = '01Ants-12'
surface_sample_num = 2000
dilation = 0.02
# inner_points = "voronoi"
inner_points = "random"
max_time_SCP = 100 # in second
```
Run
```angular2html
python Coverage_Axis_mesh.py
```
The outputs are placed in the folder `output`.
- `mesh_inner_points.obj` contains the candidate inner points.
- `mesh.obj` contains the input mesh.
- `mesh_samples_2000.obj` contains the sampled surface points that are covered.
- `mesh_selected_inner_points.obj` contains the selected inner points.

<p align="center">
<img src="./assets/fig_results_mesh.png" 
        alt="Picture" 
        width="300" 
        height="200" 
        style="display: block; margin: 0 auto" />
</p>

You may use randomly generated points inside the volume as inner candidate points by setting `inner_points = "random"
`. Notably, we already generate a sample. If you choose to produce candidates by randomly sampling inside the shape, it can be a little time consuming.
Then run
```angular2html
python Coverage_Axis_mesh.py
```

## Point Cloud Input

We use [Fast Winding Number](https://www.dgp.toronto.edu/projects/fast-winding-numbers/) for Inside-outside determination for point cloud inputs.

Please use the following commands for building the modified [libigl](https://libigl.github.io/tutorial/) at https://github.com/Frank-ZY-Dou/libigl_CA. Note that [Eigen](http://eigen.tuxfamily.org/index.php?title=Main_Page#Download) is needed for libigl; make sure you have installed it.

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

git clone https://github.com/Frank-ZY-Dou/libigl_CA.git
cd libigl_CA/
mkdir build
cd build
cmake ../
make -j8
```
Once finished, an executable file `FastWindingNumber_CA` will be generated in the folder `bin`. You can run it by
```
cd bin
./FastWindingNumber_CA ../../../input/01Ants-12_mesh.off ../../../input/01Ants-12_pc.obj ../../../input/01Ants-12_pc_random.obj
```

A point cloud will be saved to `01Ants-12_pc.obj` in the folder `input`. A randomly generated candidate skeletal points will be written to `01Ants-12_pc_random.obj` under the folder `input`. 

Then run
```angular2html
python Coverage_Axis_pc.py
```

The outputs are placed in the folder `output`.
- `pc_inner_points.obj` contains the candidate inner points.
- `pc_samples.obj` contains the points of the point cloud that is covered. **SCP is an NP-hard problem; make sure the number of to-be-covered samples is not that large.** 
- `pc_selected_inner_points.obj` contains the selected inner points.

<p align="center">
<img src="./assets/fig_results_pc.png" 
        alt="Picture" 
        width="300" 
        height="200" 
        style="display: block; margin: 0 auto" />
</p>

*Remark: We generate the point cloud inputs and inside candidates based on Fast Winding Number. The candidates are generated by randomly sampling inside the volume.
Other sampling strategies, like Voronoi-based sampling, can also be used. The core code for sampling the point cloud and generating inside candidates are given in*
```angular2html
./libigl_CA/tutorial/FastWindingNumber_CA/main.cpp
```


# More Information
### Solve Coverage Axis in MATLAB
**update:**
The code in MATLAB can be found in `./MATLAB`.

Run
```commandline
Pole_Selection_offset
```
The results are written to `./MATLAB/outputs`

You could check with output file like `vis_MA_init01Ants-27_scale0.03.obj`.



The original optimization is solved by MATLAB. In this repo, we solve SCP by Scipy in Python. *I found the solver of MILP in scipy is a little unstable compared with the MATLAB one; please suggest if you have a more powerful solver or any idea for this. Thanks ;)*
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



# Citation
```angular2html
@inproceedings{dou2022coverage,
  title={Coverage Axis: Inner Point Selection for 3D Shape Skeletonization},
  author={Dou, Zhiyang and Lin, Cheng and Xu, Rui and Yang, Lei and Xin, Shiqing and Komura, Taku and Wang, Wenping},
  booktitle={Computer Graphics Forum},
  volume={41},
  number={2},
  pages={419--432},
  year={2022},
  organization={Wiley Online Library}
}
```

# References
- https://libigl.github.io/tutorial/  Many thanks to the contributors of libigl :)
- https://www.cgal.org/
- https://gist.github.com/dendenxu/ee5008acb5607195582e7983a384e644





# fvm2d

fvm2d is a C++ program that uses a positivity-preserving finite volume method to solve the 2D radiation belt diffusion equation with a full diffusion matrix. It solves the 2D Fokker-Planck equation in ${(\alpha, \log(p))}$ coordinates and allows for the adjustment of parameters to modify the mesh, range, or time step.

With slight modification, the code can be used to solve any Fokker-Planck type equation that can be written in the conservative form.

## Install

fvm2d involves several open source C++ packages, so you need to install or add them to your environment before compiling or running the program. The packages are: 

[Eigen](https://eigen.tuxfamily.org), [xtensor](https://github.com/xtensor-stack/xtensor), [xtl](https://github.com/xtensor-stack/xtl).

One easy approach would be to simply add these libraries to the ```source``` folder of fvm2d. 

## Changes to the Makefile

If you put the above open source libraries to the ```source``` folder, then you could comment out the following line in the Makefile

```
LOCAL_INCLUDE = /Users/xtao/local/include
```
and change 
```
CCFLAGS = -Wall -Wno-class-memaccess -O2 -I$(LOCAL_INCLUDE)
```
to 
```
CCFLAGS = -Wall -Wno-class-memaccess -O2 
```

If you put the libraries in a folder whose path is PATH, then you need to modify "LOCAL_INCLUDE" accordingly. 

## Compile

After this, you may generate the executable (fvm2d) using  

```C++
make
```
If nothing happens or if you encounter any error when running the generated code, you could clean up the old files before compiling the code by using 

```C++
make clean
```

then you can run it as 

```C++
./fvm2d
```

The default input parameter file is "p.ini". 

## Introduction

fvm2d is used to solve the 2D diffusion equation in the form

$$
\frac{\partial f}{\partial t} = \sum_{i,j}\frac{1}{G}\frac{\partial}{\partial Q_i}(GD_{Q_iQ_j}\frac{\partial f}{\partial Q_j}),
$$

which can be used to describe the relativistic electron flux in Earth’s outer radiation belt. And the mathematical tool employed is Positivity-Preserving Finite Volume (PPFV) method developed by [Gao and Wu](http://epubs.siam.org/doi/10.1137/140972470), which preserves the positivity of solution and does not impose restrictions on the time step according to the CFL condition.

As for the construction of the code, the main computational and solving components are located in the **source** folder. Here we read the parameters from **p.ini** (default) and diffusion coefficients from **D** folder, construct mesh, solve the equation and finally output the files to the **output** folder. Additionally, the **plot** folder contains Python scripts that plot your results. Currently, the **cmp_ay.py** is a function that compares the default output with previous results of Albert and Young (2005) GRL results. 

## Parameter

As mentioned in part Introduction, the diffusion coefficients are stored in **D** folder, you can add your diffusion coefficients to the folder and make the appropriate changes in **p.ini** (default), or create a new ini file (for example **new.ini**) in the same path to adjust parameters such as the reference path, mesh dense, solving domain range and time step. Detailed information can be found in the defalut **p.ini**. If you choose to create a new ini file, the running command will be as follows:

```C++
./fvm2d new.ini
```

to use "new.ini" as the input parameter file.

For time dependent diffusion coefficients, boundary conditions, you will need to modify the corresponding source code.

## THINGS TO NOTE:
-- The default version of the fvm2d is to compare the fvm2d results with that of Albert and Young, GRL, 2005. The corresponding is that 

$$
f(\alpha_0 = \alpha_{0,\text{LC}}) = 0.
$$

It is also possible to change the boundary condition to 

$$
\left.\frac{\partial f}{\partial \alpha_0}\right|_{\alpha_0 = 0} = 0.
$$

In this case, the input diffusion coefficients and initial conditions of f are both changed. For initial f, see **BCs.h**. For D, you need to modify **p.ini** so that **alpha0_min_D = 0** instead of 1. Do not expect agreement for this case, due to the change to the initial f.

-- The input D's are the same as the ones used by Albert and Young (2005) GRL, and all three D's have dimension [p^2]/[t]. To use this kind of D, we processed D accordingly (see **D.cc**) before providing it to the main solver. Depending on your D, you might need to modify the **D.h** and **D.cc** class.

## Contributing to fvm2d

If you have any suggestions, please contact Peng Peng at pp140594 "AT" mail.ustc.edu.cn or Xin Tao at xtao "AT" ustc.edu.cn. 

## License

The code is open source under the MIT license. If you find it useful or base your research on it, we would appreciate it if you could cite the following paper:

Peng, P., Tao, X., Peng, Z., Jiang, Y., Gao, Z., Yang, D., et al. (2024). Modeling radiation belt dynamics using a positivity-preserving finite volume method on general meshes. Journal of Geophysical Research: Space Physics, 129, e2024JA032919. https://doi.org/10.1029/2024JA032919



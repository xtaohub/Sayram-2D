# fvm2d

fvm2d is a C++ program that uses a positivity-preserving finite volume method to solve the 2D radiation belt diffusion equation with a full diffusion matrix. It solves the 2D Fokker-Planck equation in ${(\alpha, \log(E))}$ coordinates and allows for the adjustment of parameters to modify the mesh, range, or time step.

With slight modification, the code can be used to solve any Fokker-Planck type equation that can be written in the conservative form.

## Install

fvm2d involves several open source C++ packages, so you need to install or add them to your environment before compiling or running the program. The packages are: 

[Eigen](https://eigen.tuxfamily.org), [xtensor](https://github.com/xtensor-stack/xtensor), [xtl](https://github.com/xtensor-stack/xtl), [hdf5](https://github.com/HDFGroup/hdf5)

The first three libraries are header only, so one easy approach would be to simply add them to the ```source``` folder of fvm2d.  However, you need to setup the hdf5 library before using fvm2d, because fvm2d uses HDF5 for input of D and output by default. 

## Changes to the Makefile
There is a section in the Makefile, marked by ```# ------------``` that you could change according to the setup of your coding environment.
```
# -----------------
CC = g++
LOCAL_INCLUDE = /Users/xtao/local/include
HDF5_INCLUDE = /opt/local/include
HDF5_LIB = /opt/local/lib
OPENMP =
OPT = -O2
# -----------------
```

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

The default input parameter file is "p.ini", which deals with the case in the paper by Albert and Young, GRL, 2005. 

## Introduction

fvm2d is used to solve the 2D diffusion equation in the form

$$
\frac{\partial f}{\partial t} = \sum_{i,j}\frac{1}{G}\frac{\partial}{\partial Q_i}(GD_{Q_iQ_j}\frac{\partial f}{\partial Q_j}),
$$

which can be used to describe the relativistic electron flux evolution in Earth’s outer radiation belt due to local wave particle interactions. And the mathematical tool employed is Positivity-Preserving Finite Volume (PPFV) method developed by [Gao and Wu](http://epubs.siam.org/doi/10.1137/140972470), which preserves the positivity of solution and does not impose restrictions on the time step according to the CFL condition.

As for the construction of the code, the main computational and solving components are located in the **source** folder. Here we read the parameters from **p.ini** (default) and diffusion coefficients from **D** folder, construct mesh, solve the equation and finally output the files to the **output** folder. Additionally, the **plot** folder contains Python scripts that plot your results. Currently, the **cmp_ay.py** is a function that compares the default output with previous results of Albert and Young (2005) GRL results. 

## Parameter

If you have a new configuration, called new.ini, about your case, then you can run the code as

```C++
./fvm2d new.ini
```

to use "new.ini" as the input parameter file.

For time dependent diffusion coefficients, boundary conditions, you will need to modify the corresponding source code, see the ***Equation*** class and ***Cases*** folder.

## THINGS TO NOTE:
-- The default version of the fvm2d is to compare the fvm2d results with that of Albert and Young, GRL, 2005. The corresponding is that 

$$
f(\alpha_0 = \alpha_{0,\text{LC}}) = 0.
$$

It is also possible to change the boundary condition to 

$$
\left.\frac{\partial f}{\partial \alpha_0}\right|_{\alpha_0 = 0} = 0.
$$

-- The input D's are the same as the ones used by Albert and Young (2005) GRL, and all three D's have dimension [p^2]/[t]. 

## Contributing to fvm2d

If you have any suggestions, please contact Peng Peng at pp140594 "AT" mail.ustc.edu.cn or Xin Tao at xtao "AT" ustc.edu.cn. 

## License

The code is open source under the MIT license. If you find it useful or base your research on it, we would appreciate it if you could cite the following paper:

Peng, P., Tao, X., Peng, Z., Jiang, Y., Gao, Z., Yang, D., et al. (2024). Modeling radiation belt dynamics using a positivity-preserving finite volume method on general meshes. Journal of Geophysical Research: Space Physics, 129, e2024JA032919. https://doi.org/10.1029/2024JA032919



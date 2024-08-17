# fvm2d

fvm2d is a C++ program that uses a positivity-preserving finite volume method to solve the 2D radiation belt diffusion equation with a structured mesh. It solves the 2D Fokker-Planck equation in ${(\alpha, p)}$ coordinates and allows for the adjustment of parameters to modify the mesh, range, or time step.

## Install

fvm2d involves several open source C++ packages, so you need to install or add them to your environment before compiling or running the program. The packages are: 

[Eigen](https://eigen.tuxfamily.org), [xtensor](https://github.com/xtensor-stack/xtensor), [xtl](https://github.com/xtensor-stack/xtl).

A easy way would be simply add these libraries to the ```source``` folder of fvm2d, or you could add the library to the include path of c++ compilers.

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

## Contributing to fvm2d

If you have any suggestions, please contact Peng Peng at pp140594 "AT" mail.ustc.edu.cn or Xin Tao at xtao "AT" ustc.edu.cn. 

If you base your research on this code, please consider citing the following paper:

Peng et al. (2024), Modeling radiation belt dynamics using a positivity-preserving finite volume method on general meshes, JGR: Space Physics, doi: XXXXXX

## License

Copyright (c) 2024 Xin Tao 

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the “Software”), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED “AS IS”, WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.



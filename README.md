# fvm2d

fvm2d is a C++ program that employs a positivity-preserving finite volume method to simulate the 2D diffusion of phase space density of electrons in structured mesh. It solves the 2D Fokker-Planck equation in ${(\alpha, p)}$ coordinates and allows for the adjustment of parameters to modify the mesh, range, or time step.

## Install

fvm2d involves several C++ packages, so you need to install or add them to your environment before compiling or running the program. The packages and their corresponding websites are as follows:

[Eigen](https://eigen.tuxfamily.org),

[xtensor](https://github.com/xtensor-stack/xtensor),

[xtl](https://github.com/xtensor-stack/xtl).

After the preparation, fvm2d can be simply run after the following compilation command:

```C++
make
```

then you can run it as usual C++ program by

```C++
./fvm2d
```

Notice that we don't have parameter input right now, this part will be discussed in later Parameter part.

## Introduction

fvm2d is committed to calculate the numerical results of the 2D Fokker-Planck equation

$$
\frac{\partial f}{\partial t} = \sum_{i,j}\frac{1}{G}\frac{\partial}{\partial Q_i}(GD_{Q_iQ_j}\frac{\partial f}{\partial Q_j}),
$$

which can be used to describe the relativistic electron flux in Earth’s outer radiation belt. And the mathematical tool employed is Positivity-Preserving Finite Volume (PPFV) method developed by [Gao and Wu](http://epubs.siam.org/doi/10.1137/140972470), which preserves the monotonicity and positivity of diffusion results and does not impose restrictions on the time step according to the CFL condition.

As for the construction of the code, the main computational and solving components are located in the **source** folder. Here we read the parameters from **p.ini** (default) and diffusion coefficients from **D** folder, construct mesh, solve the equation and finally output the files to the **output** folder. Additionally, the **plot** folder contains Python scripts that implement the drawing functions,performing sampling and comparing the results for energies of 0.5MeV and 2MeV at 0.1 day and 1.0 day (by default).

## Parameter

As mentioned in part Introduction, the diffusion coefficients are stored in **D** folder, you can add the relating diffusion coefficients you need to the folder and make the appropriate changes in **p.ini** (default), or create a new ini file (for example **new.ini**) in the same path to adjust parameters such as the reference path, mesh dense, solving domain range and time step. Detailed information can be found in the defalut **p.ini**. If you choose to create a new ini file, the running command will be as follows:

```C++
./fvm2d new.ini
```

to involve a second variable input to declare that a different file is chosen to be the parameter source. 

## Contributing to fvm2d

We love pull requests from everyone, and we'd be grateful that if you have any advice, bug report, and improvement ideas to share with us.

You can contact to Peng by sending an email to pp140594@mail.ustc.edu.cn

## License

TODO

# fvm2d

fvm2d is a C++ program using positivity-preserving finite volume method to simulate 2D diffusion of phase space density of electrons in structured mesh by solving 2D Fokker-Planck equation in ${(\alpha, p)}$ coordinates, which supports the input of parameters to change different mesh, range or time step.

## Install

fvm2d involves several C++ package, so it is necessary to install or add them to your environment before compiling or running the program. The packages and corresponding websites are as follows,

[Eigen](https://eigen.tuxfamily.org),

[xtensor](https://github.com/xtensor-stack/xtensor),

[xtl](https://github.com/xtensor-stack/xtl).

After the preparation, fvm2d can be simply run after compiling command:

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

which can be adopt to describe the relativistic electron flux in Earth’s outer radiation belt. And the math tool we use is Positivity-Preserving Finite Volume (PPFV) method developed by [Gao and Wu](http://epubs.siam.org/doi/10.1137/140972470), which preserve the monotone and positivity of diffusion results and the time step is not subject to the CFL condition.

As for the construction of the code, the main calculating and solving parts are in the **source** folder, in which we read the parameters from **p.ini** (default) and diffusion coefficients from **D** folder, construct mesh, solve the equation and finally output the files to **output** folder. And the **plot** folder implements the drawing function by Python, sampling and comparing the results in 0.5MeV and 2MeV at 0.1 day and 1.0 day (by default).

## Parameter

As mentioned in part Introduction, the diffusion coefficients are stored in **D** folder, you can add the relating diffusion coefficients you need to the folder and make several change in **p.ini** (default), or make a new ini file (for example **new.ini**) in the same path to change parameters including the reference path, mesh dense, solving domain range and time step. Detailed information can be referred in the defalut **p.ini**. If you choose to make a new ini file, the running command will be like:

```C++
./fvm2d new.ini
```

to involve a second variable input to declare that a different file is chosen to be the parameter source. 

## Contributing to fvm2d

We love pull requests from everyone, and we'd be grateful that if you have any advice, bug report, and improvement ideas to share with us.

You can contact to Peng by sending an email to pp140594@mail.ustc.edu.cn

## License

TODO

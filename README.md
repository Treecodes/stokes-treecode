# stokes-treecode

This repo contains files needed to run different versions of the treecode for Stokes kernels. There are two versions of the treecode and each one has two test cases. Currently all the parameters are hard-coded, so these may need to be configured before building the executables.

1. Taylor treecode for singular Stokeslet and stresslet

codes: 3D_SingStokes_Taylor_case_1.cpp, 3D_SingStokes_Taylor_case_2.cpp  
data: data_Taylor_treecode  
preprint: arxiv.org/abs/1811.12498

2. kernel-independent treecode for regularized Stokeslet

codes: KITC_MRS.cpp  
data: data.txt  
preprint: arxiv.org/abs/1902.02250

## Treecode parameters

Users set the following treecode parameters.  
P: Lagrange interpolation degree, default is 3  
numpars: number of particles, default is 10000  
N0: maximum leaf size, default is 2000  
sq_theta: square of MAC parameter, default is 0.49  
DEL: regularized Stokeslet parameter, default value is 0.02  

The following CMake macros can be defined to overide the default parameters:  
USER_PARAM_P  
USER_PARAM_NUM_PARS  
USER_PARAM_N0  
USER_PARAM_SQ_THETA  
USER_PARAM_DEL  

## Configuring and building treecode with CMake

Example configure and build commands on UWM Mortime cluster:  
```
module load icc/15.2
git clone --branch KITC_RStokeslet https://github.com/Treecodes/stokes-treecode.git
cd stokes-treecode
mkdir build
cd build
CXXFLAGS="-O2 -fopenmp" CXX=icpc cmake -D USER_PARAM_P=8 -D USER_PARAM_DEL=0.05 ..
make
```

If you don't have CMake installed, instructions to build and install one can be found at the following link:  
https://cmake.org/pipermail/cmake/2009-December/033873.html

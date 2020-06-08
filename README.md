# stokes-treecode

This repo contains files needed to run different versions of the treecode for Stokes kernels. There are two versions of the treecode and each one has two test cases. Currently all the parameters are hard-coded, so these may need to be configured before building the executables.

1. Taylor treecode for singular Stokeslet and stresslet

codes: 3D_SingStokes_Taylor_case_1.cpp, 3D_SingStokes_Taylor_case_2.cpp  
data: data_Taylor_treecode  
preprint: arxiv.org/abs/1811.12498

2. kernel-independent barycentric Lagrange treecode for regularized Stokeslets

code: KITC_MRS.cpp  
data: KITC_MRS_data.txt  
preprint: arxiv.org/abs/1902.02250

## Treecode parameters

Users set the following treecode parameters. 
n: degree of Lagrange interpolating polynomial, default is 3  
N: number of particles, default is 10000  
N0: maximum leaf size, default is 2000  
theta: MAC parameter, default is 0.7  
eps: regularized Stokeslet parameter, default is 0.02  

The following CMake macros can be defined to overide the default parameters.  
USER_PARAM_n  
USER_PARAM_N  
USER_PARAM_N0  
USER_PARAM_theta  
USER_PARAM_eps  

## Configuring and building treecode with CMake

Example configure and build commands on a Linux/Mac workstation or laptop:  
```
git clone --branch KITC_RStokeslet https://github.com/Treecodes/stokes-treecode.git
cd stokes-treecode
mkdir build
cd build
CXXFLAGS="-O2 -fopenmp" CXX=g++ cmake -D USER_PARAM_n=8 -D USER_PARAM_eps=0.05 ..
make
```

Example configure and build commands on UWM Mortimer cluster:  
```
module load icc/15.2
git clone --branch KITC_RStokeslet https://github.com/Treecodes/stokes-treecode.git
cd stokes-treecode
mkdir build
cd build
CXXFLAGS="-O2 -fopenmp" CXX=icpc cmake -D USER_PARAM_n=8 -D USER_PARAM_eps=0.05 ..
make
```

If you don't have CMake installed, instructions to build and install one can be found at the following link.  
https://cmake.org/pipermail/cmake/2009-December/033873.html  

On Mac OS X, you can install various versions of GCC via Homebrew. The default installation directory is /usr/local/Cellar/gcc.  
GCC 5: brew install gcc@5  
GCC 6: brew install gcc@6  
GCC 7: brew install gcc@7  
GCC 8: brew install gcc@8  
GCC 9: brew install gcc@9  

To use a specific version of GCC by default, say 7.2.0, first create symbolic links gcc/g++ for gcc-7/g++-7:  
```
cd /usr/local/Cellar/gcc/7.2.0/bin
ln -s gcc-7 gcc
ln -s g++-7 g++
```

Next, add the following lines to ~/.bash_profile:  
```
export PATH="/usr/local/Cellar/gcc/7.2.0/bin:$PATH"
export DYLD_LIBRARY_PATH="/usr/local/Cellar/gcc/7.2.0/lib/gcc/7:$DYLD_LIBRARY_PATH"
```

To use other versions, slightly modify above commands with corresponding version numbers.

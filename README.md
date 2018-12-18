# stokes-treecode

This repo contains files needed to run different versions of the treecode for Stokes kernels. There are two versions of the treecode and each one has two test cases. Currently all the parameters are hard-coded, so these may need to be edited before building the executables.

1. Taylor treecode for singular Stokeslet and stresslet

codes: 3D_SingStokes_Taylor_case_1.cpp, 3D_SingStokes_Taylor_case_2.cpp  
data: data_Taylor_treecode

2. kernel-independent treecode for regularized Stokeslet and rotlet

codes: 3D_RegStokes_case_1.cpp, 3D_RegStokes_case_2.cpp  
data: to be added soon

To build the executables, run the Makefile using the `make` command. Currently the Makefile is set to use the icpc compiler.


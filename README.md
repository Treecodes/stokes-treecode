# stokes-treecode

This repo contains files needed to run different versions of the treecode for Stokes kernels. There are two versions of the treecode and each one has two test cases. Currently all the parameters are hard-coded, so these may need to be edited before building the executables. To build the executables, run the Makefile using the `make` command. Currently the Makefile is set to use the icpc compiler.

This  1. barycentric Lagrange kernel-independent treecode for regularized Stokeslet and rotlet

codes: 3D_RegStokes_case_1.cpp, 3D_RegStokes_case_2.cpp  
data: to be added  
paper: L. Wang, R. Krasny, S. Tlupova (2020) A kernel-independent treecode based on barycentric Lagrange interpolation, Communications in Computational Physics 28, 1415-1436

2. Taylor treecode for singular Stokeslet and stresslet

codes: 3D_SingStokes_Taylor_case_1.cpp, 3D_SingStokes_Taylor_case_2.cpp  
data: data_Taylor_treecode   
paper: L. Wang, S. Tlupova, R. Krasny (2019) A treecode algorithm for 3D Stokeslets and stresslets, Advances in Applied Mathematics and Mechanics 11, 737-756

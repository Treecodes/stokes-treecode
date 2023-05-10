# stokes-treecode

This repo contains files needed to run different versions of the treecode for Stokes-related kernels. There are two versions of the treecode (BLTC, TTC). Currently all parameters are hard-coded. To build the executable, run the Makefile using the `make` command. Currently the Makefile is set to use the icpc compiler. The SLDM requires downloading and building LAPACK (https://netlib.org/lapack/) which used for spectral factorization of the Lanczos matrix.

1. Spectral Lanczos Decomposition Method (SLDM) with barycentric Lagrange treecode (BLTC)

**code**: SLDM_BLTC.cpp  
**data files**: ...  
**reference**: ...  

2. barycentric Lagrange treecode (BLTC) for regularized Stokeslet and rotlet

**code**: 3D_RegStokes_case_1.cpp, 3D_RegStokes_case_2.cpp  
**data files**: to be added  
**reference**: L. Wang, R. Krasny, S. Tlupova (2020) A kernel-independent treecode based on barycentric Lagrange interpolation, Communications in Computational Physics 28, 1415-1436

3. Taylor treecode (TTC) for singular Stokeslet and stresslet

**code**: 3D_SingStokes_Taylor_case_1.cpp, 3D_SingStokes_Taylor_case_2.cpp  
**data files**: data_Taylor_treecode   
**reference**: L. Wang, S. Tlupova, R. Krasny (2019) A treecode algorithm for 3D Stokeslets and stresslets, Advances in Applied Mathematics and Mechanics 11, 737-756  
**notes**: In data_Taylor_treecode there are files for two test cases,
(1) Stokeslet and stresslet particles on a sphere with system size N = 20480, 81920, 327680, 1310720,
(2) Stokeslets randomly distributed in cubes with system size N = 125K, 400K, 1000K.
Files with "lambda" in filename are weights.

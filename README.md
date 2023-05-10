# stokes-treecode

This repo contains files needed to run different versions of the treecode for Stokes-related kernels in 3D. There are two versions of the treecode, (1) BLTC (barycentric Lagrange treecode), (2) TTC (Taylor treecode). Currently all parameters are hard-coded. To build the executable, run the Makefile using the `make` command. Currently the Makefile is set to use the icpc compiler. The most recent addition uses the SLDM (Spectral Lanczos Decomposition Method) to compute correlated random displacements for Brownian dynamics simulations. The SLDM requires downloading and building LAPACK (https://netlib.org/lapack/) which is used for spectral factorization of the Lanczos matrix.

1. Spectral Lanczos Decomposition Method (SLDM) with barycentric Lagrange treecode (BLTC)

**code**: code_SLDM_BLTC.cpp  
**data files**: data_SLDM_BLTC  
**reference**: in preparation  

2. barycentric Lagrange treecode (BLTC) for regularized Stokeslet and rotlet

**code**: code_3D_RegStokes_case_1.cpp, code_3D_RegStokes_case_2.cpp  
**data files**: to be added  
**reference**: L. Wang, R. Krasny, S. Tlupova (2020) A kernel-independent treecode based on barycentric Lagrange interpolation, Communications in Computational Physics 28, 1415-1436

3. Taylor treecode (TTC) for singular Stokeslet and stresslet

**code**: code_3D_SingStokes_Taylor_case_1.cpp, code_3D_SingStokes_Taylor_case_2.cpp  
**data files**: data_Taylor_treecode   
**reference**: L. Wang, S. Tlupova, R. Krasny (2019) A treecode algorithm for 3D Stokeslets and stresslets, Advances in Applied Mathematics and Mechanics 11, 737-756  
**notes**: In data_Taylor_treecode there are files for two test cases,
(1) Stokeslet and stresslet particles on a sphere with system size N = 20480, 81920, 327680, 1310720,
(2) Stokeslets randomly distributed in cubes with system size N = 125K, 400K, 1000K.
Files with "lambda" in filename are weights.

Code used for running the Full-scale basis graphical lasso. The model is mostly implemented in MATLAB with a small amount of R code used for generating results from LatticeKrig (https://github.com/NCAR/LatticeKrig).

The run script is Run scripts/exampleScript.m. Running the model as done in the paper requires that the packages NeedMat (https://github.com/minjay/NeedMat), spherepts (https://github.com/gradywright/spherepts), and Spherical-Harmonic-Transform (https://github.com/polarch/Spherical-Harmonic-Transform) be included on the MATLAB path. These packages are required to generate the needlets, and more information on using NeedMat is contained in A Note on Spherical Needlets (http://arxiv.org/abs/1508.05406).

If you wish to use your own basis functions, these packages are unnecessary. 

Currently, a new function called fullScaleBGLXXX.m is used for each model of Z_2. This is easily changed, and I will do so as soon as I am able.

The code for logdet.m comes from 
Dahua Lin (2025). Safe computation of logarithm-determinat of large matrix (https://www.mathworks.com/matlabcentral/fileexchange/22026-safe-computation-of-logarithm-determinat-of-large-matrix), MATLAB Central File Exchange.

If you use the code, please cite the paper Modeling Large Nonstationary Spatial Data
with the Full-Scale Basis Graphical Lasso (URL coming soon
# Spline Estimation of FPCS via Manifold Conjugate Gradient Algorithm 

We propose a conjugate gradient algorithm over the product manifold to estimate FPCs. This algorithm exploits the manifold geometry structure of the overall parameter space, thus
improving its search efficiency and estimation accuracy. In addition, a roughness penalization can be easily incorporated into the algorithm with a potentially better fit. Simulation studies and the analysis of the real dataset demonstrate the appealing numerical
performance of the proposed method.

# Environment 

To reproduce the results in the paper, The user should have the following packages installed on the computer: 

* CRAN packages: snowfall, sm, mvtnorm, fpca, Rcpp, RcppArmadillo, RcppGSL
* packages built by ourselves: rOptManifold, mFPCA, FDABasics

I have uploaded the compressed files of rOptManifold, FDABasics, mFPCA. The user needs to have gsl installed before installing FDABasics. An easy way to install gsl on an OSX system is using MacPorts 
```
sudo port install gsl 
```
Caveat: part of our codes are based on the package `fpca` but it is no longer available for R 4.2.0 and higher. The "easy" and "prac" settings we used in our paper require the `easy` and `prac` data in `fpca` package. The user may download the data from CRAN.

# Reproducibility 

* PengPaul.sh, mainFPCA.R - Table 1, 2, 3, 4 

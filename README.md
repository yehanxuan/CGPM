# Spline Estimation of FPCS via Manifold Conjugate Gradient Algorithm 

We propose a conjugate gradient algorithm over the product manifold to estimate FPCs. This algorithm exploits the manifold geometry structure of the overall parameter space, thus
improving its search efficiency and estimation accuracy. A roughness penalization can be easily incorporated into the algorithm with a potentially better fit. The appealing numerical
performance of the proposed method is demonstrated by simulation
studies and the analysis of the real dataset.

# Environment 

To reproduce the results in the paper, The user should have following packages installed on the computer: 

* CRAN packages: snowfall, sm, mvtnorm, Rcpp, RcppArmadillo, RcppGSL
* packages built by ourselves: rOptManifold, mFPCA, FDABasics

I have uploaded the compressed files of rOptManifold, FDABasics, mFPCA. The user need to have gsl installed before installing FDABasics. An easy way to install gsl on OSX system is using macport 
```
sudo port install gsl 
```

# Reproducibility 

* PengPaul.sh, mainFPCA.R - Table 1, 2, 3, 4 


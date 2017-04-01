# MCMCprecision: Precision for discrete parameters in transdimensional MCMC

The R package `MCMCprecision` estimates the precision of the posterior model probabilities in transdimensional Markov chain Monte Carlo methods (e.g., reversible jump MCMC or product-space MCMC). This is useful for applications of transdimensional MCMC such as model selection, mixtures with varying numbers of components, change-point detection, capture-recapture models, phylogenetic trees, variable selection, and for discrete parameters in MCMC output in general.

To install `MCMCprecision` from GitHub, paste the following code to R (dependencies need to be installed manually):

```r
### Dependencies:
# install.packages(c("combinat", "devtools","RcppProgress","RcppArmadillo", "RcppEigen"))

library(devtools)
install_github("danheck/MCMCprecision")
```

To compile C++ code, Windows requires [Rtools](https://cran.r-project.org/bin/windows/Rtools/) and Mac [Xcode Command Line Tools](https://www.maketecheasier.com/install-command-line-tools-without-xcode/), respectively. Moreover, on Mac, it might be necessary to install the library `gfortran` manually by typing the following into the console ([required to compile the package `RcppArmadillo`](http://thecoatlessprofessor.com/programming/rcpp-rcpparmadillo-and-os-x-mavericks-lgfortran-and-lquadmath-error/)):

```
curl -O http://r.research.att.com/libs/gfortran-4.8.2-darwin13.tar.bz2
sudo tar fvxz gfortran-4.8.2-darwin13.tar.bz2 -C /
```

## Reference

* Heck, D. W., Overstall, A. M., Gronau, Q. F., & Wagenmakers, E.-J. (2017). Quantifying uncertainty in transdimensional Markov chain Monte Carlo using discrete Markov models. [arxiv:1703.10364](https://arxiv.org/abs/1703.10364)

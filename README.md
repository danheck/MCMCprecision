# MCMCprec: Precision for discrete parameters in transdimensional MCMC

The R package `MCMCprec` estimates the precision of the posterior model probabilities in transdimensional Markov chain Monte Carlo methods (e.g., reversible jump MCMC or product-space MCMC). This is useful for applications of transdimensional MCMC such as model selection, mixtures with varying numbers of components, change-point detection, capture-recapture models, phylogenetic trees, variable selection, and for discrete parameters in MCMC output in general.

To install `MCMCprec` from GitHub, paste the following code to R (dependencies need to be installed manually):

```
### Dependencies:
# install.packages(c("combinat", "sirt","devtools","RcppProgress","RcppArmadillo"))

library(devtools)
install_github("danheck/MCMCprec")
```

To compile C++ code, Windows requires [Rtools](https://cran.r-project.org/bin/windows/Rtools/) and Mac [Xcode Command Line Tools](https://www.maketecheasier.com/install-command-line-tools-without-xcode/), respectively. Moreover, on Mac, it might be necessary to install the library `gfortran` manually by typing the following into the console ([required to compile the package `RcppArmadillo`](http://thecoatlessprofessor.com/programming/rcpp-rcpparmadillo-and-os-x-mavericks-lgfortran-and-lquadmath-error/)):

```
curl -O http://r.research.att.com/libs/gfortran-4.8.2-darwin13.tar.bz2
sudo tar fvxz gfortran-4.8.2-darwin13.tar.bz2 -C /
```

## Reference

* Heck, D. W., Gronau, Q. F., Overstall, A. M., & Wagenmakers, E.-J. (2017). Estimating the Precision of Transdimensional Markov Chain Monte Carlo Methods.

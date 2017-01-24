# MCMCprec

The R package `MCMCprec` estimates the precision of the posterior model probabilities in transdimensional Markov chain Monte Carlo methods (e.g., reversible jump MCMC or product-space MCMC). This is useful for applications of transdimensional MCMC such as model selection, mixtures with varying numbers of components, change-point detection, capture-recapture models, phylogenetic trees, variable selection, and for discrete parameters in MCMC output in general.

To install `MCMCprec` from GitHub, paste the following code to R (dependencies need to be installed manually):

```
### Dependencies:
# install.packages(c("LaplacesDemon", "sirt"))
# install.packages("devtools")

library(devtools)
install_github("danheck/MCMCprec")
```

## Reference

* Heck, D. W., Gronau, Q., Overstall, A., & Wagenmakers, E.-J. (2017). Estimating the Precision of Transdimensional Markov Chain Monte Carlo Methods.

#' MCMCprec: Precision for discrete parameters in transdimensional MCMC
#'
#' MCMCprec estimates the precision of the posterior model probabilities in transdimensional Markov chain Monte Carlo methods (e.g., reversible jump MCMC or product-space MCMC). This is useful for applications of transdimensional MCMC such as model selection, mixtures with varying numbers of components, change-point detection, capture-recapture models, phylogenetic trees, variable selection, and for discrete parameters in MCMC output in general.
#'
#'
#' @author Daniel W. Heck
#' @docType package
#' @importFrom sirt dirichlet.mle
# @importFrom LaplacesDemon rdirichlet
#' @importFrom parallel clusterExport makeCluster parSapply clusterEvalQ clusterExport stopCluster
# @importFrom rARPACK eigs
#' @importFrom Matrix Matrix
#' @importFrom stats sd quantile rgamma na.omit
#' @importFrom combinat combn
#' @useDynLib MCMCprec
#' @references
#' Heck, D. W., Gronau, Q., Overstall, A. M., & Wagenmakers, E.-J. (2017).
#' Estimating the Precision of Transdimensional Markov Chain Monte Carlo Methods.
#'
"_PACKAGE"



#' Get matrix of observed transition frequencies
#'
#' Summarizes a sequence of discrete values by the observed transition frequencies.
#'
#' @param z vector of model indices (numerical or character)
#' @param labels fixed labels for models that should be included in transition matrix, e.g., \code{labels=1:20} or \code{c("m1","m2",...)}
#' @return a square matrix with transition frequencies
#' @examples
#' P <- matrix(c(.9,.1,0,
#'               .1,.6,.3,
#'               .2,.3,.5), 3, byrow=TRUE)
#' z <- sim.mc(1000, P)
#' table.mc(z, labels=1:5)
#' @export
table.mc <- function(z,
                     labels){
  if(class(z) %in% c("list", "mcmc.list")){
    z <- do.call("cbind",z)
  }
  if(missing(labels) || is.null(labels))
    labels <- sort(unique(as.vector(z)))

  if (is.matrix(z)){
    ### multiple chains
    tabs <- lapply(split(z,seq(ncol(z))),
                   table.mc, labels=labels)
    tab <- Reduce("+", tabs)
    return(tab)

  }else{
    ### single chain
    n <- length(z)
    tab <- table(factor(z[1:(n-1)], levels=labels),
                 factor(z[2:n], levels=labels))
    tab <- matrix(c(tab), nrow=length(labels),
                  dimnames=list(from=labels, to=labels))
    return(tab)
  }
}


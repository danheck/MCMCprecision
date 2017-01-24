rdirichlet <- function (n, a){
  a <- rbind(a)
  M <- ncol(a)

  if (n > dim(a)[1])
    a <- matrix(a, n, M, byrow = TRUE)

  x <- matrix(rgamma(M * n, a), ncol = M)
  sums <- x %*% rep(1, M)

  x/as.vector(sums)
}

# #' Print summary for posterior samples of stationary distribution
# #'
# #' See \code{\link{summary.stationary}}.
# #'
# #' @param ... ignored
# #' @param digits digits for rounding
# # #' @describeIn summary print summary for posterior samples


# #' @export
# print.stationary <- function(x, digits=NULL, ...){
#   if(!is.null(digits))
#     res <- lapply(summary(x), round, digits=digits)
#   else
#     res <- summary(x)
#
#   print(res,...)
# }
#
# #' @export
# head.stationary <- function(x, ...){
#   class(x) <- "matrix"
#   head(x, ...)
# }
#
# #' @export
# tail.stationary <- function(x, ...){
#   class(x) <- "matrix"
#   tail(x, ...)
# }

#' Plot restrictions
#'
#' @export
#' @param fx handle to function in question
#' @param rx handle to function that returns a sequence of feasible sequence t and x + tv
#' @param nRow number of row plots
#' @param nCol number of column plots


plot_restrictions <- function(fx, rx, nRow=3, nCol=3) {
  par(mfrow=c(nRow, nCol))
  for (i in 1:(nRow*nCol)) {
    rand_rest <- rx()
    t <- rand_rest$t
    Z <- rand_rest$Z
    nt <- length(t)
    g <- double(nt)
    for (i in 1:nt) {
      g[i] <- fx(Z[[i]])
    }
    plot(t,g,type='l')
  }
}




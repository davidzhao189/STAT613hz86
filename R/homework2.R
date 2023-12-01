#' Ridge Regression
#'
#' \code{ridge_regression} returns the ridge regression coefficient estimates
#' for a sequence of regularization parameter values.
#'
#' @param y response variables
#' @param X design matrix
#' @param lambda vector of tuning parameters
#' @export
ridge_regression <- function(y, X, lambda) {
  n <- length(y)
  B <- sapply(lambda, function(z) solve(t(X) %*% X + z*n*diag(ncol(X))) %*% t(X) %*% y)
  return(B)
}

#' Generalized Cross Validation
#'
#' \code{gcv} returns the leave-one-out
#' for a sequence of regularization parameter values.
#'
#' @param y response variables
#' @param X design matrix
#' @param lambda vector of tuning parameters
#' @export
gcv <- function(y, X, lambda) {
  diaD_square <- (svd(X)$d)^2
  n <- length(y)
  gcv <- sapply(lambda, function(z) sum(((y - X %*% solve(t(X) %*% X + n*z*diag(ncol(X))) %*% t(X) %*% y)/(1-sum(1-(n * z)/(diaD_square + n * z))/n))^2)/n)
  return(gcv)
}


#' Leave One Out
#'
#' \code{leave_one_out} returns the leave-one-out prediction error estimates
#' for a sequence of regularization parameter values.
#'
#' @param y response variables
#' @param X design matrix
#' @param lambda vector of tuning parameters
#' @export
leave_one_out <- function(y, X, lambda) {
  LOO <- rep(0,length(lambda))
  X_svd <- svd(X)
  n <- length(y)
  for (k in 1:n){
    yk_residual <- sapply(lambda, function(z) (y[k]- (X %*% solve(t(X) %*% X+z*n*diag(ncol(X))) %*% (t(X)%*%y))[k])/(1-sum((X_svd$d)^2/((X_svd$d)^2+n*z)*(X_svd$u[k,])^2)))
    LOO <- LOO + (yk_residual)^2
    }
  LOO <- LOO/n
  return(LOO)
  }

#' K-fold Cross Validation
#'
#' \code{k_fold_cv} returns K-fold cross-validation prediction error estimates
#' for a sequence of regularization parameter values.
#'
#' @param y response variables
#' @param X design matrix
#' @param k number of folds
#' @param lambda vector of tuning parameters
#' @param seed seed for random number generator
#' @export
k_fold_cv <- function(y, X, lambda, k=5, seed=12345) {
  set.seed(seed) #<- Always set a seed for reproducibility
  folds <- createFolds(y, k = k)
  
  error_sum <- rep(0,length(lambda))
  for (i in 1:k){
    y_test <- y[folds[[i]], ]
    y_train <- y[-folds[[i]], ]
    X_test <- X[folds[[i]], ]
    X_train <- X[-folds[[i]], ]
    
    beta <- ridge_regression(y_train,X_train,lambda)
    error <- apply((y_test-X_test %*% beta),2,function(z) sum(z^2))
    error_sum <- error_sum+error
  }
  return(error_sum/k)
  
}







#' Gradient Step
#'
#' @param gradf handle to function that returns gradient of objective function
#' @param x current parameter estimate
#' @param t step-size
#' @export
gradient_step <- function(gradf, t, x) {
  x_plus <- x - t*gradf(x)
  return(x_plus)
}


#' Gradient Descent (Fixed Step-Size)
#'
#' @param fx handle to function that returns objective function values
#' @param gradf handle to function that returns gradient of objective function
#' @param x0 initial parameter estimate
#' @param t step-size
#' @param max_iter maximum number of iterations
#' @param tol convergence tolerance
#' @export
gradient_descent_fixed <- function(fx, gradf, t, x0, max_iter=1e2, tol=1e-3) {
  x <- x0
  obj_value <- vector()
  rc_iterate <- vector()
  grad_value <- vector()
  rc_function <- vector()
  
  for(i in 1:max_iter){
    x_plus <- gradient_step(gradf,t,x)
    obj_value <- c(obj_value,fx(x))
    rc_iterate <- c(rc_iterate,norm(x_plus-x,type="2")/norm(x,type="2"))
    grad_value <- c(grad_value,norm(gradf(x), type="2"))
    rc_function <- c(rc_function,abs(fx(x_plus)-fx(x))/abs(fx(x)))
    if(abs(fx(x_plus)-fx(x))<= tol){
      break
    }
    x <- x_plus
  }
  
  values<-list("final_iterate_value" = x, "obj_value" = obj_value, "grad_value"=grad_value,
               "rc_function"=rc_function, "rc_iterate"=rc_iterate)
  return(values)
}

#' Backtracking
#'
#' @param fx handle to function that returns objective function values
#' @param x current parameter estimate
#' @param t current step-size
#' @param df the value of the gradient of objective function evaluated at the current x
#' @param alpha the backtracking parameter
#' @param beta the decrementing multiplier
#' @export
backtrack <- function(fx, t, x, df, alpha=0.5, beta=0.9) {
  alpha_k <- alpha
  while (fx(x-alpha_k*df(x)) >=fx(x)-0.5*alpha_k*(norm(df(x),type="2")^2)) {
    alpha_k <- beta * alpha_k
  }
  t <- alpha_k
  return(t)
}





#' Gradient Descent (Backtracking Step-Size)
#'
#' @param fx handle to function that returns objective function values
#' @param gradf handle to function that returns gradient of objective function
#' @param x0 initial parameter estimate
#' @param max_iter maximum number of iterations
#' @param tol convergence tolerance
#' @export
gradient_descent_backtrack <- function(fx, gradf, x0, max_iter=1e2, tol=1e-3) {
  x <- x0
  t <- 0.5
  obj_value <- vector()
  rc_iterate <- vector()
  grad_value <- vector()
  rc_function <- vector()
  
  for(i in 1:max_iter){
    t <- backtrack(fx, t, x, gradf, alpha=0.5, beta=0.9)
    x_plus <- gradient_step(gradf,t,x)
    obj_value <- c(obj_value,fx(x))
    rc_iterate <- c(rc_iterate,norm(x_plus-x,type="2")/norm(x,type="2"))
    grad_value <- c(grad_value,norm(gradf(x), type="2"))
    rc_function <- c(rc_function,abs(fx(x_plus)-fx(x))/abs(fx(x)))
    if(abs(fx(x_plus)-fx(x))<= tol){
      break
    }
    x <- x_plus
  }
  
  values<-list("final_iterate_value" = x, "obj_value" = obj_value, "grad_value"=grad_value,
               "rc_function"=rc_function, "rc_iterate"=rc_iterate)
  return(values)
}


#' Gradient Descent
#'
#' @param fx handle to function that returns objective function values
#' @param gradf handle to function that returns gradient of objective function
#' @param x0 initial parameter estimate
#' @param t step-size
#' @param max_iter maximum number of iterations
#' @param tol convergence tolerance
#' @export
gradient_descent <- function(fx, gradf, x0, t=NULL, max_iter=1e2, tol=1e-3) {
  if (is.null(t)==TRUE){
    values <- gradient_descent_backtrack(fx, gradf, x0, max_iter=max_iter, tol=tol)
  }else{
    values <- gradient_descent_fixed(fx, gradf, t, x0, max_iter=max_iter, tol=tol)
  }
  return(values)
}


#' Objective Function for Logistic Regression
#'
#' @param y binary response
#' @param X design matrix
#' @param beta regression coefficient vector
#' @param lambda regularization parameter
#' @export
fx_logistic <- function(y, X, beta, lambda=0) {
  L = 0
  for (i in 1:length(y)){
    L_i <- -y[i]*X[i,] %*% beta+log(1+exp(X[i,] %*% beta))
    L = L+L_i
  }
  L = L+ 0.5*lambda*(norm(beta,type="2"))^2
  return(as.numeric(L))
}



#' Gradient for Logistic Regession
#'
#' @param y binary response
#' @param X design matrix
#' @param beta regression coefficient vector
#' @param lambda regularization parameter
#' @export
gradf_logistic <- function(y, X, beta, lambda=0) {
  G = rep(0,length(beta))
  for (i in 1:length(y)){
    G_i <- as.numeric((-y[i] + log(1+exp(X[i,] %*% beta))/(1+log(1+exp(X[i,] %*% beta))))) * X[i,]
    G = G+G_i
  }
  G = G+ lambda*beta
  return(G)
}





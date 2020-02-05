update_L <- function(L, Y, X, rho, tol, maxit, quiet)
{
  p <- nrow(L)
  k <- ncol(L)
  
  # Elements above diagonal are zero, leading elements positive
  positive_idx <- 1
  if(k >= 2){
    positive_idx <- c(positive_idx, 1 + cumsum(p:(p - k + 2)))  
  }
  
  
  obj_f <- function(x){
    L2 <- matrix(0, nrow = p, ncol = k)
    L2[lower.tri(L2, diag = TRUE)] <- x
    return(obj_rcpp(L2, Y, X, rho))
  }
  
  grad_f <- function(x){
    L2 <- matrix(0, nrow = p, ncol = k)
    L2[lower.tri(L2, diag = TRUE)] <- x
    L2 <- t(jac_rcpp(L2, Y, X, rho))
    return(as.vector(L2[lower.tri(L2, diag = TRUE)]))
  }
  

  lower <- rep(-Inf, p * k - (k - 1) * k / 2)
  upper <- rep(Inf, p * k -  (k - 1) * k / 2)
  lower[positive_idx] <- 0
  out <- stats::optim(par = as.vector(L[lower.tri(L, diag = TRUE)]), fn = obj_f,
                      gr = grad_f, method = "L-BFGS-B", lower = lower, upper = upper,
                      control = list(trace = 1 * !quiet, maxit = maxit,
                                     factr = tol / .Machine$double.eps))
  if(out$convergence != 0){
    warning("L-BFGS-B did not converge")
  }
  L[upper.tri(L)] <- 0
  L[lower.tri(L, diag = TRUE)] <- out$par
  return(L)
}
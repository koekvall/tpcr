update_L <- function(L, Y, X, rho, tol, maxit, quiet)
{
  p <- nrow(L)
  k <- ncol(L)
  obj_f <- function(x){
    L2 <- matrix(x, nrow = p, ncol = k)
    return(obj_rcpp(L2, Y, X, rho))
  }
  
  grad_f <- function(x){
    L2 <- matrix(x, nrow = p, ncol = k)
    return(as.vector(t(jac_rcpp(L2, Y, X, rho))))
  }
  
  out <- stats::optim(par = as.vector(L), fn = obj_f, gr = grad_f, method = "BFGS",
               control = list(trace = 1 * !quiet, maxit = maxit, abstol = tol,
                              reltol = tol))
  if(out$convergence != 0){
    warning("BFGS did not converge")
  }
  return(matrix(out$par, p, k))
}
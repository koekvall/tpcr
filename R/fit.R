jlpcr <- function(Y, X, k, rho = 0, tol = 1e-10, maxit = 1e3, center = FALSE, quiet = TRUE, L)
{
  # Define constants
  p <- ncol(X)
  n <- nrow(X)
  r <- ncol(Y)
  
  if(center){
    X <- scale(X, scale = F)
    Y <- scale(Y, scale = F)
  }
  
  if(k >= min(n, p)) stop("JLPCR requires k < min(n, p)")
  if(missing(L)){
    # Use probabilistic principal components estimates
    e_S <- eigen(crossprod(X) / n, symmetric = TRUE)
    
    U <- e_S$vectors[, 1:k, drop = F]
    tau <- mean(e_S$values[(k + 1):p])
    D <- e_S$values[1:k] / tau - 1
    L <- chol_spsd(U %*% diag(D, k) %*% t(U) ,tol = 1e-8)$L
    if(ncol(L) > k) warning("PPCA starting calue has more than k columns; dropping.")
    L <- L[, 1:k, drop = F]
    rm(e_S, U, D, tau)
  }
  L <- matrix(L, ncol = k)
  ############################################################################
  # Find Estimates
  ############################################################################
  L <- update_L(L = L, Y = Y, X = X, rho = rho, tol = tol, maxit = maxit,
                  quiet = quiet)
  Z <- X %*% L
  alpha <- solve(crossprod(Z) + rho * crossprod(L), crossprod(Z, Y))
  Sigma <- crossprod(Y - Z %*% alpha) / n
  tau <- sum(X^2) -  sum(t(X) * (L %*% qr.solve(diag(1, k) +
                                               crossprod(L), t(L)) %*% t(X)))
  tau <- tau / (n * p)

  beta <-  L %*% alpha
  
  return(list(beta = beta,
              Sigma = Sigma, Sigma_X = tau * (diag(1, p) + tcrossprod(L)),
              L = L, tau = tau))
}

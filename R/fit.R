#' Fit principal components regression using a jointly normal likelihood
#'
#' @param Y n x r matrix of responses
#' @param X n x p matrix of predictors
#' @param k Number of principal components to use
#' @param rho Ridge penalty on alpha in representation beta = L alpha
#' @param tol Tolerance for L-BFGS-B on profile likelihood
#' @param maxit Maximum number of iterations of L-BFGS-B algorithm
#' @param center If TRUE, responses and predictors are centered by their sample mean
#' @param scale If TRUE, divide each column in X by its sample standard deviation
#' @param quiet If FALSE, print information from L-BFGS-B algorithm
#' @param L Starting value in L-BFGS-B for the Cholesky root
#'          in the decomposition Sigma_X = tau (I_p + LL^T)
#' @return List with estimates:
#'          beta (regression coefficient)
#'          Sigma (response covariance matrix)
#'          Sigma_X (predictor covariance matrix)
#'          L (Cholesky root in decomposition of Sigma_X)
#'          tau (smallest eigenvalue of Sigma_X; has multiplicity p - k)
#' @useDynLib jlpcr, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#' @importFrom Rcpp evalCpp
#' @export
jlpcr <- function(Y, X, k, rho = 0, tol = 1e-10, maxit = 1e3, center = TRUE,
  scale = FALSE, quiet = TRUE, L)
{
  # Define constants
  p <- ncol(X)
  n <- nrow(X)
  r <- ncol(Y)

  if(n < p & rho == 0){
    warning("The objective function is unbounded when n < p and rho = 0; do not expect reliable estimates!")
  }

  if(center){
    X <- scale(X, scale = F)
    Y <- scale(Y, scale = F)
  }

  if(k >= min(n, p)) stop("jlpcr requires k < min(n, p)")
  if(missing(L)){
    # Use probabilistic principal components estimates
    e_S <- eigen(crossprod(X) / n, symmetric = TRUE)

    U <- e_S$vectors[, 1:k, drop = F]
    tau <- mean(e_S$values[(k + 1):p])
    D <- e_S$values[1:k] / tau - 1
    L <- chol_spsd(U %*% diag(D, k) %*% t(U) ,tol = 1e-8)$L
    if(ncol(L) > k) warning("PPCA starting value has more than k columns; dropping.")
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

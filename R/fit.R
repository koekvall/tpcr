#' Targeted Principal Components Regression
#' 
#' \code{tpcr} fits a principal components regression by maximizing a joint
#'  multivariate normal pseud-likelihood.
#' 
#' This is the only function exported from the package with the same name. 
#' The likelihood maximized is that for n independent observations 
#' from a normal multivariate response linear regression model where the column space
#' of the p-by-r coefficient matrix is spanned by the k leading eigenvectors (corresponding to
#' the largest eigenvalues) of the predictors' covariance matrix.
#' This covariance matrix is assumed to be spiked, meaning its
#' smallest eigenvalue has multiplicity p - k.
#'
#' @param Y Matrix (n x r) of responses; rows correspond to observations
#' @param X Matrix (n x p) of predictors; rows correspond to observations
#' @param k Integer number of principal components to use; can be a vector
#' @param rho Numeric ridge penalty on alpha in representation beta = L alpha
#' @param tol Numeric tolerance for L-BFGS-B on profile log-likelihood
#' @param maxit Integer maximum number of iterations of L-BFGS-B algorithm
#' @param center_Y If TRUE, responses are centered by their sample average.
#' @param center_X If TRUE, predictors are centered by their sample average.     
#' @param scale_Y If TRUE, responses are scaled to have unit sample standard deviation.
#' @param scale_X If TRUE, predictors are scaled to have unit sample standard deviation.
#' @param quiet If TRUE, suppresses information from optim (L-BFGS-B)
#' @param L Matrix (p x k) starting value in L-BFGS-B for the Cholesky root
#'          in the decomposition Sigma_X = tau (I + LL')
#' @param m Numeric penalty in information criterion
#'          - 2 * log-likelihood + m * n_params; can be a vector
#' @param covmat If TRUE, calculates asymptotic covariance matrix of vec(beta)
#' @param Xnew Matrix of new observations to predict the response for
#' @return If length(k) = 1, returns a list with estimates and information criterion.
#'         If scale_X or scale_Y are TRUE, estimates are rescaled to original scale.
#'          
#'         If length(k) > 1, returns a list of lists of length k + 1 where for
#'         j in 1:k the jth element is the list returned by tpcr with k = k[j]
#'         and the (k + 1)th element is a vector of the same length as m where
#'         the jth element is the k selected by the IC corresponding to m[j].
#' @useDynLib tpcr, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#' @importFrom Rcpp evalCpp
#' @export
tpcr <- function(Y, X, k, rho = 0, tol = 1e-10, maxit = 1e3, center_Y = TRUE,
  center_X = TRUE, scale_Y = FALSE, scale_X = FALSE, quiet = TRUE, L, m = 2,
  covmat = FALSE, Xnew = X)
{
  #############################################################################
  # Argument checking and starting values
  #############################################################################
  stopifnot(is.matrix(X))
  k <- sort(k)
  num_k <- length(k)
  stopifnot(is.numeric(k), is.null(dim(k)), all(k == floor(k)))
  p <- ncol(X)
  n <- nrow(X)
  if(k[num_k] >= min(n, p)) stop("tpcr requires k < min(n, p)")
  
  stopifnot(is.numeric(Y))
  if(!is.matrix(Y)){
    Y <- matrix(Y, nrow = n)
  }
  r <- ncol(Y)
  stopifnot(is.numeric(rho), length(rho) == 1, rho >= 0)
  if(n < p & rho == 0){
    warning("The objective function is unbounded when n < p and rho = 0;
            do not expect reliable estimates!")
  }
  stopifnot(is.numeric(tol), length(tol) == 1, tol > 0)
  stopifnot(is.numeric(maxit), length(maxit) == 1, maxit > 0, maxit == floor(maxit))
  stopifnot(is.logical(center_Y), length(center_Y) == 1)
  stopifnot(is.logical(center_X), length(center_X) == 1)
  stopifnot(is.logical(scale_Y), length(scale_Y) == 1)
  stopifnot(is.logical(scale_X), length(scale_X) == 1)
  stopifnot(is.logical(quiet), length(quiet) == 1)
  # Get starting value
  if(missing(L) | num_k > 1){
    e_S <- eigen(crossprod(X) / n, symmetric = TRUE)
  }
  if(missing(L)){
    # Use probabilistic principal components estimates
    U <- e_S$vectors[, 1:k[1], drop = F]
    tau <- mean(e_S$values[(k[1] + 1):p])
    D <- e_S$values[1:k[1]] / tau - 1
    L <- chol_spsd(U %*% diag(D, k[1]) %*% t(U) , tol = 1e-8)$L
    L <- L[, 1:k[1], drop = F]
    rm(D, tau, U)
  } else{
    stopifnot(is.matrix(L), ncol(L) == k, nrow(L) == p)
  }
  
  stopifnot(is.numeric(m), is.null(dim(m)))
  stopifnot(is.logical(covmat), length(covmat) == 1)
  stopifnot(is.matrix(Xnew), ncol(Xnew) == p)

  #############################################################################
  # Scale data if required
  #############################################################################
  ybar <- rep(0, r)
  xbar <- rep(0, p)
  if(center_Y){
    ybar <- colMeans(Y)
    Y <- sweep(Y, 2, ybar)
  }
  if(center_X){
   xbar <- colMeans(X)
   X <- sweep(X, 2, xbar)
  }
  sy <- rep(1, r)
  if(scale_Y){
    sy <- apply(Y, 2, stats::sd)
    Y <- sweep(Y, 2, 1 / sy, FUN = "*")
  }
  sx <- rep(1, p)
  if(scale_X){
    sx <- apply(X,  2, stats::sd)
    X <- sweep(X, 2, 1 / sx, FUN = "*")
  }
  
  #############################################################################
  # Fit model for each k_i in k
  #############################################################################
  out_list <- list()
  for(ii in 1:num_k){
    # Get estimates for possibly scaled data
    L <- update_L(L = L, Y = Y, X = X, rho = rho, tol = tol, maxit = maxit,
                  quiet = quiet)
    sL <- svd(L)
    U <- sL$u
    d <- sL$d^2
    Z <- X %*% U
    gam <- solve(crossprod(Z) + diag(rho, k[ii]), crossprod(Z, Y))
    Sigma <- crossprod(Y - Z %*% gam) / n
    tau <- sum(X^2) -  sum(t(X) * (L %*% qr.solve(diag(1, k[ii]) +
                                                  crossprod(L), t(L)) %*% t(X)))
    tau <- tau / (n * p)
    d <- d * tau
    Psi <- U %*% diag(d, k[ii]) %*% t(U)
    beta <-  U %*% gam
    Sigma_X <- diag(tau, p) + Psi
    
    # Get information criteria for possibly scaled data
    n_param <- r * (r + 1) / 2 # Sigma
    n_param <- n_param + k[ii] * r # gamma
    n_param <- n_param + 1 # tau
    n_param <- n_param + k[ii] + k[ii] * p - k[ii] * (k[ii] + 1) / 2 #Psi
    if(center_Y) n_param <- n_param + r # Intercept
    if(center_X) n_param <- n_param + p # Predictor mean
    # k eigenvalues, k orthonormal vectors requires (p - 1) + (p - 2) ... (p - k)
    # = k + kp - sum_i^k i = k + kp - k(k + 1) / 2
    IC <-  sum(mvtnorm::dmvnorm(Y -  X %*% beta, sigma = Sigma, log = T))
    IC <- IC + sum(mvtnorm::dmvnorm(X, sigma = Sigma_X, log = T))
    IC <- -2 * IC + m * n_param
    
    # Get get estimates on original scale for returning
    if(scale_Y){
      Sigma <- sweep(Sigma, 2, sy, FUN = "*")
      Sigma <- sweep(Sigma, 1, sy, FUN = "*")
      beta <-  sweep(beta, 2, sy, FUN = "*")
    }
    if(scale_X){
      Sigma_X <- sweep(Sigma_X, 2, sx, FUN = "*")
      Sigma_X <- sweep(Sigma_X, 1, sx, FUN = "*")
      beta <-    sweep(beta, 1, 1 / sx, FUN = "*")
    }
    # Covariance of estimates on original scale
    C <- NULL
    if(covmat){
      if(scale_X){
        U <-  sweep(U, 1, sx, FUN = "*")
      } 
      C <- kronecker(Sigma, U %*% ((1 / (d + tau)) * t(U)))
    }
    # Prediction or fitted values on original scale
    Yhat <- sweep(sweep(Xnew, 2, xbar) %*% beta, 2, ybar, FUN = "+")
    # Return list for k_i
    out_list[[ii]] <- list(ybar = ybar, xbar = xbar, b = beta,
                           Sigma = Sigma,  Sigma_X = Sigma_X,
                           IC = IC, cov_b = C, Yhat = Yhat, L = L)
    # Starting value for next iter; use directions of first sample evecs not in span
    if(ii < num_k){
      u_new <- e_S$vectors[, 1:(k[ii + 1] - k[ii]), drop = F] - 
        sL$u %*% crossprod(sL$u, e_S$vectors[, 1:(k[ii + 1] - k[ii]), drop = F])
      u_new <- scale(u_new, center = F, scale = apply(u_new, 2,
                                                      function(x)sqrt(sum(x^2))))
      U_hat <- cbind(sL$u, u_new)
      L <- chol_spsd(U_hat %*% diag(c(sL$d^2, rep(sL$d[k[ii]]^2,
                       k[ii + 1] - k[ii])), k[ii + 1]) %*%
                       t(U_hat), tol = 1e-8)$L[, 1:k[ii + 1]]
    }
  }
  #############################################################################
  # Prepare return list
  #############################################################################
  if(num_k == 1){
    out_list <- out_list[[1]]
  } else{
    all_IC <- do.call(rbind, lapply(out_list, function(x)x$IC))
    out_list[[num_k + 1]] <- k[as.vector(apply(all_IC, 2, which.min))]
    names(out_list) <- c(paste0("fit_k_", k), "k_star")
  }
  return(out_list)
}

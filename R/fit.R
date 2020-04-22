#' Joint Likelihood Principal Components Regression
#' 
#' \code{jlpcr} fits a principal components regression by maximizing the joint likelihood 
#'  of multivariate normal responses and predictors.
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
#' @param mean_Y Bool indicating the mean of Y should be a parameter; if FALSE,
#'               it is assumed E(Y) = 0.
#' @param mean_X Bool indicating the mean of X should be a parameter; if FALSE,
#'               it is assumed E(X) = 0.           
#' @param scale Bool indicating each predictor is divided by its sample
#'              standard deviation
#' @param quiet Bool indicating no information printed from L-BFGS-B algorithm
#' @param L Matrix (p x k) starting value in L-BFGS-B for the Cholesky root
#'          in the decomposition Sigma_X = tau (I + LL')
#' @param m Numeric penalty in information criterion
#'          - 2 * log-likelihood + m * n_params; can be a vector
#' @return If length(k) = 1, returns a list with estimates:
#'          beta (regression coefficient),
#'          Sigma (response covariance matrix),
#'          Sigma_X (predictor covariance matrix),
#'          L (Cholesky root in decomposition Sigma_X = tau[I + LL']),
#'          tau (smallest eigenvalue of Sigma_X; has multiplicity p - k), and
#'          IC (information criterion).
#'          
#'         If length(k) > 1, returns a list of lists of length k + 1 where for
#'         j in 1:k the jth element is the list returned by jlpcr with k = k[j]
#'         and the (k + 1)th element is a vector of the same length as m where
#'         the jth element is the k selected by the IC corresponding to m[j].
#' @useDynLib jlpcr, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#' @importFrom Rcpp evalCpp
#' @export
jlpcr <- function(Y, X, k, rho = 0, tol = 1e-10, maxit = 1e3, mean_Y = TRUE,
  mean_X = TRUE, scale = FALSE, quiet = TRUE, L, m = 2)
{
  Y <- as.matrix(Y)
  X <- as.matrix(X)
  # Define constants
  k <- sort(k)
  p <- ncol(X)
  n <- nrow(X)
  r <- ncol(Y)

  if(n < p & rho == 0){
    warning("The objective function is unbounded when n < p and rho = 0; do not expect reliable estimates!")
  }
  
  mu_Y <- rep(0, ncol(Y))
  mu_X <- rep(0, ncol(X))
  if(mean_Y){
    Y <- scale(Y, scale = F)
    mu_Y <- attr(Y, "scaled:center")
  }
  if(mean_X){
    X <- scale(X, scale = F)
    mu_X <- attr(X, "scaled:center")
  }
  
  if(scale) X <- scale(X, center = F)
  
  if(max(k) >= min(n, p)) stop("jlpcr requires k < min(n, p)")
  if(missing(L)){
    # Use probabilistic principal components estimates
    e_S <- eigen(crossprod(X) / n, symmetric = TRUE)

    U <- e_S$vectors[, 1:k[1], drop = F]
    tau <- mean(e_S$values[(k[1] + 1):p])
    D <- e_S$values[1:k[1]] / tau - 1
    L <- chol_spsd(U %*% diag(D, k[1]) %*% t(U) , tol = 1e-8)$L
    L <- L[, 1:k[1], drop = F]
    rm(D, tau, U)
  }
  L <- matrix(L, ncol = k[1])
  
  ############################################################################
  # Cycle through k
  ############################################################################
  if(length(k) > 1){
    if(missing(e_S)){
      e_S <- eigen(crossprod(X) / n, symmetric = TRUE)
    }
    out_list <- list()
    ii <- 1
    for(k_i in k){
      out_list[[ii]] <- jlpcr(Y = Y, X = X, k = k_i, rho = rho, tol = tol,
                              maxit = maxit, mean_Y = mean_Y, mean_X = mean_X,
                              scale = scale, quiet = quiet, L = L, m = m)
      
      # New starting value; use directions of first sample evecs not in span
      if(ii == length(k)) break #  No starting values for next iter
      s_L <- svd(out_list[[ii]]$L)
      u_new <- e_S$vectors[, 1:(k[ii + 1] - k_i), drop = F] - 
        s_L$u %*% crossprod(s_L$u, e_S$vectors[, 1:(k[ii + 1] - k_i), drop = F])
      u_new <- scale(u_new, center = F, scale = apply(u_new, 2, function(x)sqrt(sum(x^2))))
      U_hat <- cbind(s_L$u, u_new)
      L <- chol_spsd(U_hat %*% diag(c(s_L$d^2, rep(s_L$d[k_i]^2, k[ii + 1] - k_i)), k[ii + 1]) %*%
                       t(U_hat), tol = 1e-8)$L[, 1:k[ii + 1]]
      ii <- ii + 1
    }
    all_IC <- do.call(rbind, lapply(out_list, function(x)x$IC))
    out_list[[ii + 1]] <- k[as.vector(apply(all_IC, 2, which.min))]
    names(out_list) <- c(paste0("fit_k_", k), "k_star")
    return(out_list)
  }
  
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
  Sigma_X <- tau * (diag(1, p) + tcrossprod(L))
  
  n_param <- r * (r + 1) / 2 # Cond covmat
  n_param <- n_param + k * r # alpha
  n_param <- n_param + 1 # tau
  n_param <- n_param + k + k * p - k * (k + 1) / 2 #SPSD rank k
  if(mean_Y) n_param <- n_param + r # Intercept
  if(mean_X) n_param <- n_param + p # Predictor mean
  # k eigenvalues, k orthonormal vectors requires (p - 1) + (p - 2) ... (p - k)
  # = k + kp - sum_i^k i = k + kp - k(k + 1) / 2
  IC <-  sum(mvtnorm::dmvnorm(Y -  X %*% beta, sigma = Sigma, log = T))
  IC <- IC + sum(mvtnorm::dmvnorm(X, sigma = Sigma_X, log = T))
  IC <- -2 * IC + m * n_param
  
  return(list(mu_Y = mu_Y, mu_X = mu_X, beta = beta,
              Sigma = Sigma, Sigma_X = tau * (diag(1, p) + tcrossprod(L)),
              L = L, tau = tau, IC = IC))
}

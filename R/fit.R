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
#'          in the decomposition Sigma_X = tau (I + LL')
#' @param m Penalty in information criterion
#'          - 2 * log-likelihood + m * n_params; can be a vector
#' @return List with estimates:
#'          beta (regression coefficient),
#'          Sigma (response covariance matrix),
#'          Sigma_X (predictor covariance matrix),
#'          L (Cholesky root in decomposition Sigma_X = tau[I + LL']),
#'          tau (smallest eigenvalue of Sigma_X; has multiplicity p - k), and
#'          IC (information criterion)
#' @useDynLib jlpcr, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#' @importFrom Rcpp evalCpp
#' @export
jlpcr <- function(Y, X, k, rho = 0, tol = 1e-10, maxit = 1e3, center = TRUE,
  scale = FALSE, quiet = TRUE, L, m = 2)
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
  
  if(k[1] >= min(n, p)) stop("jlpcr requires k < min(n, p)")
  if(missing(L)){
    # Use probabilistic principal components estimates
    e_S <- eigen(crossprod(X) / n, symmetric = TRUE)

    U <- e_S$vectors[, 1:k[1], drop = F]
    tau <- mean(e_S$values[(k[1] + 1):p])
    D <- e_S$values[1:k[1]] / tau - 1
    L <- chol_spsd(U %*% diag(D, k[1]) %*% t(U) , tol = 1e-8)$L
    if(ncol(L) > k[1]) warning("PPCA starting value has more than k columns; dropping.")
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
      out_list[[ii]] <- jlpcr(Y = Y, X = X, k = k_i, rho = rho, tol = tol, maxit = maxit,
            center = center, scale = scale, quiet = quiet, L = L, m = m)
 
      # New starting value; use directions of first sample evecs not in span
      if(ii == length(k)) break
      s_L <- svd(out_list[[ii]]$L)
      u_new <- e_S$vectors[, 1:(k[ii + 1] - k_i), drop = F] - 
        s_L$u %*% crossprod(s_L$u, e_S$vectors[, 1:(k[ii + 1] - k_i), drop = F])
      u_new <- scale(u_new, center = F, scale = apply(u_new, 2, function(x)sqrt(sum(x^2))))
      U_hat <- cbind(s_L$u, u_new)
      L <- chol_spsd(U_hat %*% diag(c(s_L$d^2, rep(s_L$d[k_i]^2, k[ii + 1] - k_i)), k[ii + 1]) %*%
                       t(U_hat), tol = 1e-8)$L[, 1:k[ii + 1]]
      ii <- ii + 1
    }
    out_list[[ii + 1]] <- k[which.min(unlist(lapply(out_list, function(x)x$IC[1])))]
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
  n_param <- n_param + 2 * p * k - k^2 - (p - k - 1) 
  
  IC<-  sum(mvtnorm::dmvnorm(Y -  X %*% beta, sigma = Sigma, log = T))
  IC <- IC + sum(mvtnorm::dmvnorm(X, sigma = Sigma_X, log = T))
  IC <- -2 * IC + m * n_param
  
  return(list(beta = beta,
              Sigma = Sigma, Sigma_X = tau * (diag(1, p) + tcrossprod(L)),
              L = L, tau = tau, IC = IC))
}

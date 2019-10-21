set.seed(4)

# Dimensions (k < p & k < n)
n <- 200; p <- 10; k <- 2; r <- 1

# Error variances
ssy0 <- 1; tau0 <- 1

Sigma0 <- diag(ssy0, r)

L0 <- matrix(0, p, k)
L0[, 1] <- runif(p)
L0[2:p, 2] <- runif(p - 1)

SigmaX0 <- tau0 * (diag(p) + tcrossprod(L0))

# True coefficients
alpha0 <- matrix(runif(r * k, min = -1, max = 1), k, r)
beta0 <- L0 %*% alpha0


# Data
X <- matrix(rnorm(n * p), nrow = n, ncol = p) %*% chol(SigmaX0)
Y <- X %*% beta0 + matrix(rnorm(n * r), n, r) %*% chol(Sigma0)

# PCR estimate
e_X <- eigen(crossprod(X) / n)
tau_pcr <- mean(e_X$values[(k + 1):p])
L_pcr <- e_X$vectors[, 1:k] %*% diag(sqrt(e_X$values[1:k]/tau_pcr - 1), k)
Sigma_X_pcr <- tau_pcr * (diag(p) + tcrossprod(L_pcr))

# Our estimate
test_obj <- function(x)
{
  L <- matrix(x, p, k)
  return(mlpcr:::obj_rcpp(L, Y, X, 1))
}

grad_f <- function(x)
{
  numDeriv::grad(test_obj, x)
}

grad_my <- function(x){
  L <- matrix(x, p, k)
  return(c(t(mlpcr:::jac_rcpp(L, Y, X, 1))))
}

# 
# grad_f(out$par)
# 
# L_hat <- matrix(0, p, k)
# L_hat[, 1] <- out$par[1:p]
# L_hat[2:p, 2] <- out$par[-c(1:p)]
# 
# tau_hat <- sum(diag(X %*% qr.solve(diag(p) + tcrossprod(L_hat), t(X)))) / (n * p)
# 
# Sigma_X_hat <- tau_hat * (diag(p) + tcrossprod(L_hat))
# 
# alpha_hat <- solve(crossprod(X %*% L_hat), crossprod(X %*% L_hat, Y))
# beta_hat <- L_hat %*% alpha_hat
# Sigma_hat <- crossprod(Y - X %*% beta_hat) / n
# 
# res_mat[ii, ] <- c(norm(Sigma_X_hat - SigmaX0) , norm(Sigma_X_pcr - SigmaX0), norm(crossprod(X)/n - SigmaX0))
# }

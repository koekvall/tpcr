set.seed(4)

# Dimensions (k < p & k < n)
n <- 50; p <- 40; k <- 5; r <- 2

# Error variances
ssy0 <- 1; tau0 <- 1

Sigma0 <- diag(ssy0, r)

L0 <- matrix(rnorm(p * k), p, k)

SigmaX0 <- tau0 * (diag(p) + tcrossprod(L0) / 5)

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

# A small simulation
n_sims <- 50
res_mat <- matrix(0, n_sims, 4)
for(ii in 1:n_sims){
  X <- matrix(rnorm(n * p), nrow = n, ncol = p) %*% chol(SigmaX0)
  Y <- X %*% beta0 + matrix(rnorm(n * r), n, r) %*% chol(Sigma0)
 
  beta_mlpcr <- mlpcr:::mlpcr(Y, X, k, rho = 0)$beta
  beta_pcr <- coef(pls::pcr(Y ~ 0 + X, ncomp = k))[, , 1]
  beta_pls <- coef(pls::plsr(Y ~ 0 + X, ncomp = k, method = "oscores"))[, , 1]
  beta_env <- Renvlp::xenv(X, Y, u = k, asy = FALSE)$beta
  
  res_mat[ii, ] <- c(norm(beta_mlpcr - beta0, "2"), norm(beta_pcr - beta0, "2"),
                     norm(beta_pls - beta0, "2"),  norm(beta_env - beta0, "2"))
}


# Advertising
ad_dat <- read.csv("~/Desktop/Advertising.csv")
Y <- matrix(ad_dat$sales[1:100], ncol = 1)
X <- as.matrix(ad_dat[1:100, 2:4])
Y_test <- matrix(ad_dat$sales[101:200], ncol = 1)
X_test <- as.matrix(ad_dat[101:200, 2:4])
beta_mlpcr <- mlpcr:::mlpcr(Y, X, 1, rho = 0)$beta
beta_unscaled <- mlpcr:::mlpcr(Y, X, 1, stand = F)$beta

beta_pcr <- coef(pls::pcr(Y ~ X, ncomp = 1))[, , 1]
beta_pls <- coef(pls::plsr(Y ~ X, ncomp = 1, method = "oscores"))[, , 1]
beta_env <- Renvlp::xenv(X, Y, u = 1, asy = FALSE)$beta
beta_ols <- coef(lm(Y ~ X))

rmse_mlpcr <- norm(Y_test - X_test %*% beta_mlpcr, "F")
rmse_pcr <- norm(Y_test - X_test %*% beta_pcr, "F")
rmse_pls <- norm(Y_test - X_test %*% beta_pls, "F")
rmse_env <- norm(Y_test - X_test %*% beta_env, "F")
rmse_ols <- norm(Y_test - X_test %*% beta_ols, "F")


# FREDMD
setwd("~/GDrive/Research/ts_pls/code/")
library(tidyverse)
library(pls)
library(tseries)
library(doParallel)
library(forecast)

fred <- read.csv("../FRED/050119.csv")
fred[, 1] <- as.character(fred[, 1])
transforms <- as.numeric(read_csv("../FRED/fred_transforms.csv"))

# Transform codes in Appendix of https://s3.amazonaws.com/real.stlouisfed.org/wp/2015/2015-012.pdf
transform_fred <- function(x, type){
  if(type == 1) return(x)
  else if(type == 2) return(c(NA, diff(x)))
  else if(type == 3) return(c(NA, NA, diff(diff(x))))
  else if(type == 4) return(log(x))
  else if(type == 5) return(c(NA, diff(log(x))))
  else if(type == 6) return(c(NA, NA, diff(diff(log(x)))))
  else if(type == 7) return(c(NA, NA, diff(x[-1] / x[-length(x)] - 1)))
  else(stop("type must be in 1:7"))
}

# Transform data
for(ii in 2:129){
  fred[, ii] <- transform_fred(fred[, ii], transforms[ii - 1])
}

# Remove series with many missing values
fred <- as_tibble(fred) %>% select(-ACOGNO, -ANDENOx, -OILPRICEx, -TWEXMMTH, -UMCSENTx)
fred <- fred[complete.cases(fred), ]
train_idx <- 1:floor(0.8 * nrow(fred))
train_dat <- fred[train_idx, ]

# Standardize 
train_dat[, -1] <- scale(train_dat[, -1])
n_train <- nrow(train_dat)



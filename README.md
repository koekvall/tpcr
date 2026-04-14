# tpcr (Targeted Principal Components Regression)

An R package for multivariate principal components regression. Instead of
the usual two-step procedure (estimate principal components, then regress),
`tpcr` jointly estimates the principal components and regression coefficients
by maximizing a normal pseudo-likelihood. This targets the components most
relevant for prediction.

## Installation

```r
install.packages("devtools")
devtools::install_github("koekvall/tpcr")
```

## Example

```r
library(tpcr)

# Simulate data: n observations, p predictors, r responses
n <- 200
p <- 50
r <- 3
k <- 2  # number of principal components

# Spiked covariance: two strong components plus noise
U <- svd(matrix(rnorm(p * k), p, k))$u
Sigma_X <- U %*% diag(c(10, 5)) %*% t(U) + diag(1, p)
X <- MASS::mvrnorm(n, rep(0, p), Sigma_X)

# Coefficient matrix with columns in the span of the leading eigenvectors
gamma <- matrix(rnorm(k * r), k, r)
beta <- U %*% gamma
Y <- X %*% beta + matrix(rnorm(n * r), n, r)

# Fit targeted PCR with k = 2
fit <- tpcr(Y, X, k = 2)

# Estimated coefficients
fit$b

# Fitted values
head(fit$Yhat)

# Select k by information criterion (AIC and BIC)
fit_select <- tpcr(Y, X, k = 1:5, m = c(2, log(n)))
fit_select$k_star  # selected k for each penalty
```

## Main function

The package exports a single function, `tpcr(Y, X, k, ...)`, with arguments:

- `Y` -- response matrix (n x r)
- `X` -- predictor matrix (n x p)
- `k` -- number of principal components (integer or vector for model selection)
- `rho` -- ridge penalty (default 0)
- `m` -- penalty for information criterion (default 2 for AIC; use `log(n)` for BIC)
- `covmat` -- if `TRUE`, returns the asymptotic covariance matrix of the estimated coefficients
- `Xnew` -- new predictor matrix for prediction

See `?tpcr` for the full list of arguments.

## Citation

The corresponding paper is published in the [Journal of Multivariate Analysis](https://www.sciencedirect.com/science/article/pii/S0047259X22000318) (also available on [arXiv](https://arxiv.org/abs/2004.14009)).

A BibTeX entry is

```bibtex
@article{ekvall2022targeted,
  title={Targeted principal components regression},
  author={Ekvall, Karl Oskar},
  journal={Journal of Multivariate Analysis},
  volume={190},
  pages={104973},
  year={2022},
  publisher={Elsevier},
  doi={10.1016/j.jmva.2022.104973}
}
```

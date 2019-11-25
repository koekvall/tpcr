chol_spsd <- function(A, tol = 1e-15)
{
  m <- ncol(A)
  L <- matrix(0, m, m)
  L[1, 1] <- sqrt(A[1, 1])
  r <- 1
  counter <- 0
  d <- 1
  for(ii in 2:m){
    for(jj in 1:(ii - counter - 1)){
      L[ii, jj] <- (A[ii, d[jj]] - sum(L[d[jj], 1:(jj - 1)] * L[ii, 1:(jj - 1)])) / L[d[jj], jj]
    }
    v <- sqrt(abs(A[ii, ii] - sum(L[ii, 1:(ii - counter - 1)]^2)))
    if(v^2 <= tol){
      counter <- counter + 1
    }else{
      L[ii, ii - counter] <- v
      d <- c(d, ii)
      r <- r + 1
    }
  }
  return(list(L = L[, 1:r, drop = FALSE], r = r, d = d))
}
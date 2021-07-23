#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
double obj_rcpp(arma::mat L, arma::mat Y, arma::mat X, double rho)
{
  const unsigned int n = Y.n_rows;
  const unsigned int p = X.n_cols;
  const unsigned int k = L.n_cols;

  arma::mat Z = X * L;
  arma::mat H = Z.t() * Z + rho * L.t() * L;
  
  arma::mat a = arma::solve(H, Z.t() *  Y, arma::solve_opts::likely_sympd);
  Y -= Z * a;
  double obj;
  double sign;
  arma::log_det(obj, sign, Y.t() * Y);

  arma::mat M = L * L.t();
  M.diag() += 1.0;
  double detval;
  arma::log_det(detval, sign, M);
  obj += detval;
  
  //arma::mat K = X * arma::solve(M, X.t());
  //obj += p * std::log(arma::accu(K.diag()));
  obj += p * std::log(arma::accu(X.t() % arma::solve(M, X.t())));
  obj += rho * arma::accu(arma::square(L * a));
  return obj;
}

// [[Rcpp::export]]
arma::mat jac_rcpp(arma::mat L, arma::mat Y, arma::mat X, double rho)
{
  const unsigned int n = Y.n_rows;
  const unsigned int p = X.n_cols;
  const unsigned int k = L.n_cols;
  
  arma::mat Z = X * L;
  arma::mat H = Z.t() * Z + rho * L.t() * L;
  
  //Replace H by its inverse
  arma::vec d;
  arma::eig_sym(d, H, H);
  H = H * arma::diagmat(1.0 / d) * H.t();
  
  arma::mat a = H * Z.t() *  Y;
  
  arma::mat E = Y - Z * a;

  // First term contribution; log |(Y - Za)'(Y - Za)|
  // (Y'Y)^{-1}Y' = VD^{-2}V' V D U' = VD^{-1}U'
  arma::mat U1;
  arma::mat V1;
  arma::vec s1;
  arma::svd_econ(U1, s1, V1, E);
  
  arma::mat J = - 2.0 * a * V1 * arma::diagmat(1.0 / s1) * U1.t() * X;

  // Jacobian for \bar{\alpha}; save and add up from different terms
  arma::mat A = -2.0 * V1 * arma::diagmat(1.0 / s1) * U1.t() * Z;

  // // Second term constribution; log |I_p + LL'|
  // // Callculate inverse of (I + LL')
  arma::mat U2;
  arma::mat V2;
  arma::vec s2;
  arma::svd_econ(U2, s2, V2, L);
  s2 = s2 % s2; // Replace by eigenvalues of LL'
  arma::mat Inv = -U2 * arma::diagmat(s2 / (1.0 + s2)) * U2.t();
  Inv.diag() += 1.0;
  J += 2.0 * L.t() * Inv;

  // Third term contribution; p x logtr(X'X(I + LL')^{-1})
  arma:: mat F = X * Inv;
  double c = -2.0 * p / arma::accu(X % F);
  J += c * L.t() * F.t() * F;

  // // Fourth term contribution; rho ||L a||_F^2
  J += 2.0 * rho * a * a.t() * L.t();
  A += 2 * rho * a.t() * L.t() * L;
  // 
  // From first and fourth term, d \bar{\alpha}
  
  J -= H * A.t() * Y.t() * Z * H * (Z.t() * X + L.t() * rho);
  
  J -= H * Z.t() * Y * A * H * (Z.t() * X + L.t() * rho);
  
  J += H * A.t() * Y.t() * X;
  
  return J;
}

// // [[Rcpp::export]]
// arma::mat chol_lr(arma::mat A, double tol)
// {
//   // Low rank Cholesky factorization of Canto et al. (2015)
//   const int m = A.n_cols;
//   arma::mat L(m, m, arma::fill::zeros);
//   arma::ivec c(m, arma::fill::zeros);
//   int c_length = 1;
//   int count = 0;
//   int r = 1;
//   L(0, 0) = std::sqrt(A(0, 0));
//   
//   for(int ii = 1; ii < m; ii++){
//     int idx = ii - count - 1;
//     for(int jj = 0; jj <= idx; jj++){
//       L(ii, jj) = A(ii, c(jj));
//       if(jj > 0){
//         arma::mat l1 = L.submat(c(jj), 0, c(jj), jj - 1);
//         arma::mat l2 = L.submat(ii, 0, ii, jj - 1);
//         L(ii, jj) -= arma::accu(l1 % l2);
//       }
//       L(ii, jj) *= (1.0 / L(c(jj), jj));
//     }
//   double v = std::sqrt(A(ii, ii) - arma::accu(
//               arma::square(L.submat(ii, 0, ii, idx))));
//     if(abs(v) <= tol){
//       count += 1;
//     } else{
//       L(ii, ii - count) = v;
//       c(c_length) = ii;
//       c_length += 1;
//       r += 1;
//     }
//   }
//   return L.cols(0, r - 1);
// }

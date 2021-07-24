#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
// [[Rcpp::export]]
arma::mat MLE_BL_cpp(int r, double eps, arma::mat Mx, arma::mat Bk, arma::mat M2E){

  int n = M2E.n_rows;
  double bias = 1.0;

  arma::mat Bkk(n,r);
  arma::mat Mekk(n,n);
  arma::mat EFF(r,r);
  arma::mat EZF(n,r);
  arma::mat Mxk(n,n);
  arma::mat Mek(n,n);
  arma::mat invMek(n,n);
  arma::mat invMxkBk(n,r);
  arma::mat thetakk(n,r+1);
  arma::mat D(n,1);
  arma::mat DD(r,r);
  arma::mat DR(r,1);
  D.fill(1.0);
  DR.fill(1.0);
  DD.fill(0.0);
  Mek.fill(0.0);
  invMek.fill(0.0);
  Mekk.fill(0.0);
  DD.diag() = DR;

  while(bias > eps){
    Mek.diag() = M2E;
    invMek.diag() = D/M2E;
    Mxk = Bk*Bk.t() + Mek;
    invMxkBk = invMek*Bk*pinv(DD + Bk.t()*invMek*Bk);
    EFF = invMxkBk.t()*Mx*invMxkBk + DD - Bk.t()*invMxkBk;
    EZF = Mx*invMxkBk;
    Bkk = EZF*pinv(EFF);
    Mekk.diag() = diagvec(Mx - Bkk*invMxkBk.t()*Mx);
    bias = norm(Bkk - Bk,"fro") + norm(Mekk.diag() - M2E,"fro");
    Bk = Bkk;
    M2E = Mekk.diag();
  }
  thetakk.cols(0,r-1) = Bk;
  thetakk.col(r) = M2E;
  return thetakk;
}

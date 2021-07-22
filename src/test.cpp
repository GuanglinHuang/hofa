#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;

// [[Rcpp::export]]
arma::mat MLE_BL_cpp(double UU, double HH, arma::mat YY, arma::mat FF, arma::mat ZZ){
  double h = HH;
  double u = UU;
  double pi = 3.14159265358979323846;

  arma::mat X = FF;
  int t = X.n_rows;
  arma::mat Y = YY.rows(1,t-1);
  arma::mat W(t-1,t-1);
  arma::mat XX(t-1,4);
  arma::mat uu(t-1,1);
  arma::mat ee(t-1,1);
  arma::mat hh(t-1,1);
  arma::mat two(t-1,1);
  arma::mat pipi(t-1,1);

  W.fill(0);uu.fill(u); ee.fill(1.0);hh.fill(h);two.fill(2.0);pipi.fill(pi);
  XX.col(0) = ee; XX.col(1) = X.rows(1,t-1);
  XX.col(2) = ZZ.rows(0,t-2) - uu;
  XX.col(3) = X.rows(1,t-1)%(ZZ.rows(0,t-2) - uu);
  arma::mat f = ZZ.rows(0,t-2);

  W.diag() = ((ee/sqrt(two%pipi))%exp(-((uu - f)/h)%((uu - f)/h)/two))/h;

  arma::mat Beta = (pinv(XX.t()*W*XX)*XX.t()*W*Y).t();

  return Beta;
}

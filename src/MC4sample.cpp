#include <Rcpp.h>

using namespace Rcpp;
// [[Rcpp::export]]
SEXP  MC4smaple(SEXP XX, SEXP MM2, SEXP NN, SEXP TT){
  using namespace Rcpp ;
  NumericVector X(XX);
  NumericVector M2(MM2);
  int T = as<int>(TT);
  int N = as<int>(NN);
  NumericMatrix MC4(N, pow(N,3));

  for(int i=0; i<N; i++) {
    int iiN = i * T;
    for(int j=0; j<N; j++) {
      int jjN = j * T;
      for(int k=0; k<N; k++) {
        int kkN = k * T;
        for(int l=0; l<N; l++) {
          int llN = l * T;
          // compute 4th cumulant element
          double elem = 0.0;
          for (int tt = 0; tt < T; tt++) {
            elem += X[iiN + tt]*X[jjN + tt]*X[kkN + tt]*X[llN + tt]/T;
          }
          MC4(l,i*pow(N,2)+j*N+k) = elem - M2[i+j*N]*M2[k+l*N]- M2[i+k*N]*M2[j+l*N] - M2[i+l*N]*M2[j+k*N];
        } // loop l
      } // loop k
    } // loop j
  } // loop i
  return MC4;
}

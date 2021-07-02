#' Estimating the latent factors and factor loadings in high dimensional factor model using principal component methods based on the covariance or correlation matrix.
#'
#' @param X A matrix or data frame with t rows (samples) and n columns (variables).
#' @param C Characteristics, used in Projected PCA.
#' @param r The number of factors.
#' @param scale logical. If \code{TRUE}, the variance of columns of X are normalized to 1 before factor estimation.
#' @param method Method to use: "\code{PCA}", Basic Principal Component Analysis; "\code{P-PCA}", Fan et al.(2016)'s Projected PCA.
#' @param ... Any other parameters.
#' @return The number of non-Gaussian and Gaussian factors determined by selected approach.
#' \itemize{
#'   \item{\code{f}}{  Estimated factors.}
#'   \item{\code{u}}{  Estimated factor loadings.}
#'   \item{\code{e}}{  Estimated errors.}
#'   \item{\code{m2e}}{  Diagonal elements of the covariance matrix of errors, only provided in \code{PCA} and \code{ML}.}
#'   \item{\code{rho}}{  Auto regression coefficients of errors, only provided in \code{ML-GLS}, \code{ML-ITE} and \code{ML-EM}.}
#' }
#' @examples
#' n = 100;t = 200;k = 2;
#' par_f = c(0.5,0,2,2,Inf,Inf);
#' par_e = c(1,0,2,Inf,0,0,0);
#' rho_ar = c(0.5,0.2);
#' data = hofa.sim(n,t,k,par_f,par_e,rho_ar)$X;
#' M2.mle(data,r = 2,method = "ML-EM");

n = 100
t = 10
g11 = function(x){x^3-2*x}
g12 = function(x){x^2*sin(x)}
g21 = function(x){x^2-1}
g22 = function(x){x}
g31 = function(x){x^3-x^2-x+1}
g32 = function(x){exp(-x^2)*x^2}
C = matrix(rnorm(n*2),n,2)
W = matrix(NA,n,3)
W[,1] <- g11(C[,1])+g12(C[,2])
W[,2] <- g21(C[,1])+g22(C[,2])
W[,3] <- g31(C[,1])+g32(C[,2])
FF = matrix(rnorm(t*3),t,3)
EE = matrix(rnorm(t*n),t,n)
X = W%*%t(FF) + t(EE)

TraceRatio(F_pca,FF)
TraceRatio(F_hat,FF)

TraceRatio(W_pca,W)
TraceRatio(W_pp,W)

plot(C[,1],W_pca[,1])
plot(C[,1],W_pca[,2])
plot(C[,1],W_pca[,3])

plot(C[,1],W[,1])
plot(C[,1],W[,2])
plot(C[,1],W[,3])

M2.pca = function(X,C,r,scale = F,method = c("PCA","P-PCA"),...){

  n = NCOL(X)
  t = NROW(X)
  X1 = scale(X,scale = scale)
  Mx = cov(X1)

  #PCA
  u_pca = eigen(Mx)$vectors[,1:r]*sqrt(n)
  f_pca = X1%*%u_pca/n
  c_pca = f_pca%*%t(u_pca)
  e_pca = X1 - c_pca
  if(method == "PCA"){
    con = list(f = f_pca,u = u_pca,e = e_pca)
  }
  if(method == "P-PCA"){
    X1 = scale(t(X),scale = scale)
    d = NCOL(C)
    for (i in 1:t) {
      at = paste0("s(C[,",1,"],J)")
      if(d > 1){
        for (j in 2:d){
          b = paste0("s(C[,",j,"],J)")
          at = paste0(at,"+",b)
        }
      }
      fo = paste0("X1[,i]~",at)
      reg_fit = gam::gam(formula(fo))
      PX[,i] <- reg_fit$fitted.values
    }
    eig = eigen(t(PX)%*%PX/n)
    F_hat = eig$vectors[,1:3]*sqrt(t)
    G_hat = PX%*%F_hat/t
    W_hat = X%*%F_hat/t
    Gamma_hat = W_hat-G_hat

    con = list(f = F_hat,u = W_hat,Gx = G_hat,gamma = Gamma_hat, ev = eig$values)
  }

   return(con)
}



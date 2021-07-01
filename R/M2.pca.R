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

n = 5000
t = 50
g1 = function(x){x}
g2 = function(x){x^2-1}
g3 = function(x){x^3-2*x}
C = matrix(rnorm(n*3),n,3)
W = matrix(NA,n,3)
W[,1] <- g1(C[,1])
W[,2] <- g2(C[,2])
W[,3] <- g3(C[,3])
FF = matrix(rnorm(t*3),t,3)
EE = matrix(rnorm(t*n),t,n)
X = W%*%t(FF) + t(EE)
X = scale(X,scale = F)
PX = X
for (i in 1:t) {
  reg_fit = gam(X[,i] ~ s(C[,1])+s(C[,2])+s(C[,3]))
  PX[,i] <- reg_fit$fitted.values
}

eig = eigen(t(PX)%*%PX/n)

F_hat = eig$vectors[,1:3]*sqrt(t)

GX = PX%*%F_hat/t

plot(C[,1],GX[,3])
plot(C[,2],GX[,2])
plot(C[,3],GX[,1])

eigpca = eigen(t(X)%*%X/n)
F_pca = eigpca$vectors[,1:3]*sqrt(t)
W_pca = X%*%F_pca/t

TraceRatio(F_pca,FF)
TraceRatio(F_hat,FF)

TraceRatio(GX,W)
TraceRatio(W_pca,W)

plot(C[,1],W_pca[,3])
plot(C[,2],W_pca[,2])
plot(C[,3],W_pca[,1])
plot(C[,3],W[,3])

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
  m2e_pca = diag(diag(cov(e_pca)),n,n)
  m2_pca = cov(c_pca) + m2e_pca

}



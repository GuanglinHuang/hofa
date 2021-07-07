#' Estimating the latent factors using principal component methods
#' @description Estimating the latent factors and factor loadings in high dimensional factor model using principal component methods based on the covariance or correlation matrix.
#' @param X A matrix or data frame with t rows (samples) and n columns (variables).
#' @param C Characteristics, a matrix with n rows (variables) and d columns (characteristics), used in Projected PCA.
#' @param r The number of factors.
#' @param center logical. If \code{TRUE}, the mean of columns of X are normalized to 0 before factor estimation.
#' @param scale logical. If \code{TRUE}, the variance of columns of X are normalized to 1 before factor estimation.
#' @param method Method to use: "\code{PCA}", Basic Principal Component Analysis; "\code{P-PCA}", Fan et al.(2016)'s Projected PCA.
#' @param J The number of sieve terms in Projected PCA. Default to use the criterion in Fan et al.(2016).
#' @param ... Any other parameters.
#' @return A list of factors, factor loadings and other information, see below.
#' \itemize{
#'   \item{\code{f}}{  Estimated factors.}
#'   \item{\code{u}}{  Estimated factor loadings.}
#'   \item{\code{e}}{  Estimated errors.}
#'   \item{\code{ev}}{ Eigenvalues of covariance matrix.}
#'   \item{\code{G}}{  Estimated non-parametric functions, only provided in \code{P-PCA}.}
#'   \item{\code{gamma}}{ Errors in factor loading matrix, only provided in \code{P-PCA}.}
#' }
#' @examples
#' n = 100;t = 10;d = 1;r = 3;
#' g1 = function(x){x^3-2*x};
#' g2 = function(x){x^2-1};
#' g3 = function(x){x};
#' C = matrix(rnorm(n*d),n,d);W = matrix(NA,n,r);
#' W[,1] <- g1(C);W[,2] <- g2(C);W[,3] <- g3(C);
#' FF = matrix(rnorm(t*r),t,r);
#' EE = matrix(rnorm(t*n),t,n);
#' X = W%*%t(FF) + t(EE);
#' M2.pca(t(X),C = C,r,method = "P-PCA",J = 4);




M2.pca = function(X,C = NULL,r,center = F,scale = F,method = c("PCA","P-PCA"),J = NULL,...){
  s <- gam::s
  if(method == "P-PCA" && is.null(C)){
    stop("Please input characteristic variables for Projected PCA.")
  }
  n = NCOL(X)
  t = NROW(X)

  X1 = scale(X,scale = scale)
  Mx = cov(X1)

  #PCA
  u_pca = eigen(Mx)$vectors[,1:r]*sqrt(n)
  ev = eigen(Mx)$values
  f_pca = X1%*%u_pca/n
  c_pca = f_pca%*%t(u_pca)
  e_pca = X1 - c_pca
  if(method == "PCA"){
    con = list(f = f_pca,u = u_pca,e = e_pca, ev = eigen(Mx)$values/t)
  }
  if(method == "P-PCA"){
    if(is.null(J)){
      J = round(3*(n*min(n,t))^(1/4))-1
    }
    X1 = scale(t(X),center = center,scale = scale)
    d = NCOL(C)
    PX = X1
    for (i in 1:t) {
      at = paste0("s(C[,",1,"],J)")
      if(d > 1){
        for (j in 2:d){
          b = paste0("s(C[,",j,"],J)")
          at = paste0(at,"+",b)
        }
      }
      fo = paste0("X1[,i]~",at)
      reg_fit = gam::gam(stats::formula(fo))
      PX[,i] <- reg_fit$fitted.values
    }
    eig = eigen(t(PX)%*%PX/n)
    F_hat = eig$vectors[,1:r]*sqrt(t)
    G_hat = PX%*%F_hat/t
    W_hat = X1%*%F_hat/t
    e_hat = t(X1 - W_hat%*%t(F_hat))
    Gamma_hat = W_hat-G_hat

    con = list(f = F_hat,u = W_hat,e = e_hat, ev = eig$values, G = G_hat,gamma = Gamma_hat)
  }

   return(con)
}



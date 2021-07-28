#' Estimating the latent factors using fourth-order generalized moment methods
#' @description Estimating the latent factors and factor loadings in high dimensional factor model using generalized moment methods based on the third-order cumulant.
#' @param X A matrix or data frame with t rows (samples) and n columns (variables).
#' @param r The number of factors.
#' @param kappa An integer. The weight between \code{M} and \code{VWV}, \code{M} denotes the covariance matrix and \code{VWV} denotes the GMM matrix. \code{kappa} default to 0.
#' @param initial The method used to initialize the factor loadings and the moments of error. \code{PCA} denotes Principal Component Analysis and \code{MLE} denotes Maximum Likelihood Estimation in Bai and Li(2012).
#' @param skewness Logical. If \code{TRUE}, the skewness moment structure equations will be used.
#' @param W_diag Logical. If \code{TRUE}, the weight matrix \code{W} only has nonzero diagonal elements, default to \code{FALSE}.
#' @param identity Logical. If \code{TRUE}, the moments of errors are unified to identity, default to \code{FALSE}.
#' @param delta An integer. The hard threshold value of \code{W} matrix in order to guarantee non-singularity of \code{W}, default to 1/log(t).
#' @param eps The iteration error, default to 10^-6. Available for initializing the estimators by Maximum Likelihood method.
#' @param ... Any other parameters.
#' @return A list of factors, factor loadings and other information, see below.
#' \itemize{
#'   \item{\code{f}}{  Estimated factors.}
#'   \item{\code{u}}{  Estimated factor loadings.}
#'   \item{\code{e}}{  Estimated errors.}
#'   \item{\code{ev}}{ Eigenvalues of covariance matrix.}
#' }
#' @examples
#' n = 100
#' t = 200
#' k = 2
#' par_f = c(0.5,0,2,2,Inf,Inf)
#' par_e = c(1,0,2,Inf,0,0,0)
#' rho_ar = c(0.5,0.2)
#' data = hofa.sim(n,t,k,par_f,par_e,rho_ar)$X
#' M4.gmm(data,r = 2,kappa = 0,initial = "PCA")


M4.gmm <- function(X,r,kappa = 0,initial = c("PCA","MLE"),skewness = FALSE,W_diag = FALSE,identity = FALSE,delta = NULL,eps = 10^-6,...){

  n = NCOL(X)
  t = NROW(X)

  if(skewness == T){
    m = n + n + n + 1
  }else{
    m = n + n + 1
    }

  X1 = scale(X,scale = FALSE)
  Mx = cov(X1)

  #PCA
  d_pca = eigen(Mx)$values[1:r]/n
  u_pca = eigen(Mx)$vectors[,1:r]%*%diag(sqrt(d_pca),r,r)*sqrt(n)
  f_pca = X1%*%u_pca%*%diag(1/sqrt(d_pca),r,r)/n
  c_pca = f_pca%*%t(u_pca)
  e_pca = X1 - c_pca

  if(initial == "PCA"){
    u_inl_id = u_pca
    sigma_id = colMeans(e_pca^2)
    sk_id = colMeans(e_pca^3)
    kt_id = colMeans(e_pca^4) - 3*sigma_id^2
  }
  #use MLE to initialize sigma_e
  if(initial == "MLE"){
    B_k = u_pca #initialize
    m2e =  as.matrix(diag(cov(e_pca)*(t-1)/t))

    theta_k = hofa::MLE_BL_cpp(r = r,eps = sqrt(eps),Mx = Mx,Bk = B_k,M2E = m2e)

    u_inl =  theta_k[,1:r]
    Me_inl = diag(theta_k[,r+1])
    f_inl = t(solve(t(u_inl)%*%diag(1/diag(Me_inl))%*%u_inl)%*%t(u_inl)%*%diag(1/diag(Me_inl))%*%t(X1))
    c_inl = f_inl%*%t(u_inl)
    e_inl = X1 - c_inl

    ev_inl = eigen(t(c_inl)%*%c_inl/t)

    u_inl_id = ev_inl$vectors[,1:r]
    sigma_id = diag(Me_inl)
    sk_id = colMeans(e_inl^3)
    kt_id = colMeans(e_inl^4)
  }

  sigma_e = sigma_id
  sk_e = sk_id
  kt_e = kt_id

  if(identity == T){
    sigma_e = rep(mean(sigma_e),n)
    sk_e = rep(mean(sk_e),n)
    kt_e = rep(mean(kt_e),n)
  }

  #g(x) functions
  #m1
  v1 = colMeans(X)
  #m2
  v2 = t(X1)%*%X1/t - diag(sigma_e)
  #m3
  v3 = t(X1*X1)%*%X1/t - diag(sk_e)
  #m4
  v4 = t(X1*X1*X1)%*%(X1)/t - 3*(t(X1)%*%X1/t)*(colMeans(X1^2)%*%t(rep(1,n))) - diag(kt_e)

  if(skewness == T){
    V = cbind(v1,v2,v3,v4)
  }else{
    V = cbind(v1,v2,v4)
  }

  #
  W_sol = matrix(0,m,m)
  for (tt in 1:t) {
    v1_i = X[tt,]
    v2_i = X1[tt,]%*%t(X1[tt,]) - diag(sigma_e)
    v3_i = (X1[tt,]*X1[tt,])%*%(t(X1[tt,])) - diag(sk_e)
    v4_i = (X1[tt,]*X1[tt,]*X1[tt,])%*%t(X1[tt,]) - 3*(X1[tt,]%*%t(X1[tt,]))*(colMeans(X1^2)%*%t(rep(1,n))) - diag(kt_e)
    if(skewness == T){
      V_i = cbind(v1_i,v2_i,v3_i,v4_i)
    }else{
      V_i = cbind(v1_i,v2_i,v4_i)
    }
    W_sol = W_sol + t(V_i)%*%(diag(n) - u_inl_id%*%t(u_inl_id))%*%V_i
  }
  W_sol = W_sol/t

  if(W_diag == T){
    W_sol = diag(diag(W_sol),m,m)
  }

  eig_w = eigen(W_sol)
  ev_w = as.numeric(eig_w$values[1:(min(n,t))])
  if(is.null(delta)){
    delta = 1/log(t)
  }
  id = which(ev_w > delta)
  evec = matrix(as.numeric(eig_w$vectors[,1:(min(n,t))]),m,min(n,t))
  WW = (evec[,id])%*%diag(1/ev_w[id])%*%t(evec[,id])

  eig_gmm =  eigen(kappa*t(X)%*%X/t + V%*%WW%*%t(V))

  # plot(eig_gmm$values)

  uv_gmm = as.numeric(eig_gmm$vectors[,1:r])
  u_gmm = matrix(uv_gmm,n,r)
  ev_gmm = (as.numeric(eig_gmm$values))[1:(min(n,t))]

  f_gmm = X%*%u_gmm

  e_gmm = X - f_gmm%*%t(u_gmm)

  con = list(f = f_gmm,u = u_gmm,e = e_gmm, ev = ev_gmm)

  return(con)
}


# rep = 100
# con= matrix(NA,rep,4)
# set.seed(1231)
# for (iii in 1:rep){
#   n = 50;t = 50;d = 2;r = 2;u = 2;
#   g1 = function(x){x};
#   g2 = function(x){x};
#   C = matrix(rnorm(n*d),n,d);W = matrix(NA,n,r);
#   W[,1] <- g1(C[,1]);W[,2] <- g2(C[,2]);
#   FF = matrix(NA,t,r);
#   FF[,1] = rnorm(t) - u;FF[,2] = rnorm(t) + u
#   sige = diag(sqrt(runif(n,1,5)))
#   EE = scale((matrix(rchisq(t*n,1),t,n) - 1),scale = F)%*%sige;
#   X = FF%*%t(W) + EE;
#
#   pca2 = M2.pca(X,r = 2,method = "PCA")
#   gmm2 = M2.gmm(X,r = 2,kappa = 0.5,identity = F,initial = "PCA")
#   gmm3 = M3.gmm(X,r = 2,kappa = 0.5,identity = F,initial = "PCA")
#   gmm4 = M4.gmm(X,r = 2,kappa = 0.5,identity = F,skewness = F,initial = "PCA")
#
#   # TraceRatio(gmm2$u,W)
#   # TraceRatio(gmm3$u,W)
#   # TraceRatio(gmm4$u,W)
#   # TraceRatio(pca2$u,W)
#
#  con[iii,1] = TraceRatio(gmm2$u,W)
#  con[iii,2] = TraceRatio(gmm3$u,W)
#  con[iii,3] = TraceRatio(gmm4$u,W)
#  con[iii,4] = TraceRatio(pca2$u,W)
#  print(iii)
# }
#
# boxplot(con[,c(1,2,3,4)])
# boxplot(con[,4:6])
#
# pca = M2.pca(X,r=3,method = "PCA")
# ppca = M2.pca(X,C = C,r = 3,method = "P-PCA",J = 4)
# TraceRatio(gmm$u,W)
# TraceRatio(pca$u,W)
# TraceRatio(ppca$u,W)
#
# TraceRatio(gmm$f,FF)
# TraceRatio(pca$f,FF)
# TraceRatio(ppca$f,FF)
#
# kap_tol = seq(0,5,0.1)
# rep = 100
# hgl = 1
# con_tol = vector()
# for (kappa in kap_tol){
#
#   con = vector()
#   for (iii in 1:rep) {
#     n = 30;t = 100;d = 1;r = 3;
#     g1 = function(x){x};
#     g2 = function(x){x};
#     g3 = function(x){x};
#     C = matrix(rnorm(n*d),n,d);W = matrix(NA,n,r);
#     W[,1] <- g1(C);W[,2] <- g2(C);W[,3] <- g3(C);
#     FF = (matrix(rchisq(t*r,1),t,r)-1)%*%diag(sqrt(c(3,1,0.5)));
#     EE = matrix(rchisq(t*n,1),t,n)-1;
#     MU = rep(1,t)%*%t(runif(n,0,4))
#     X = FF%*%t(W) + EE;
#
#     gmm = M2.gmm(X,r = 3,kappa = kappa,sigma = rep(1,n),initial = "PCA")
#
#     con[iii] <- as.numeric(TraceRatio(gmm$u,W))
#   }
#
#   con_tol = cbind(con_tol,con)
#
#   print(hgl)
#   hgl = hgl +1
# }
#
# plot(colMeans(con_tol),type="b")
# boxplot(con_tol)

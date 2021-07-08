#' Estimating the latent factors using generalized moment methods
#' @description Estimating the latent factors and factor loadings in high dimensional factor model using generalized moment methods based on the covariance matrix.
#' @param X A matrix or data frame with t rows (samples) and n columns (variables).
#' @param r The number of factors.
#' @param kappa An integer. The weight between \code{M} and \code{VWV}, \code{M} denotes the covariance matrix and \code{VWV} denotes the GMM matrix. \code{kappa} default to 0.
#' @param sigma A \code{n x 1} vector, the variance of the error terms. If it is \code{NULL}, the variance vector will be initialized by PCA or MLE.
#' @param initial The method used to initialize the factor loadings and the variance of errors. \code{PCA} denotes Principal Component Analysis and \code{MLE} denotes Maximum Likelihood Estimation in Bai and Li(2012).
#' @param delta An integer. The hard threshold value of \code{W} matrix in order to guarantee non-singularity of \code{W}, default to 1/log(t).
#' @param eps The iteration error, default to 10^-8. Available for initializing the estimators by Maximum Likelihood method.
#' @param ... Any other parameters.
#' @return A list of factors, factor loadings and other information, see below.
#' \itemize{
#'   \item{\code{f}}{  Estimated factors.}
#'   \item{\code{u}}{  Estimated factor loadings.}
#'   \item{\code{e}}{  Estimated errors.}
#'   \item{\code{ev}}{ Eigenvalues of covariance matrix.}
#' }
#' @examples
#' n = 100;t = 200;k = 2;
#' par_f = c(0.5,0,2,2,Inf,Inf);
#' par_e = c(1,0,2,Inf,0,0,0);
#' rho_ar = c(0.5,0.2);
#' data = hofa.sim(n,t,k,par_f,par_e,rho_ar)$X;
#' M2.gmm(data,r = 2,kappa = 0,sigma = rep(1,n),initial = "PCA");


n = 10;t = 300;d = 1;r = 3;
g1 = function(x){x^3-2*x};
g2 = function(x){x^2-1};
g3 = function(x){x};
C = matrix(rnorm(n*d),n,d);W = matrix(NA,n,r);
W[,1] <- g1(C);W[,2] <- g2(C);W[,3] <- g3(C);
FF = matrix(rnorm(t*r),t,r);
EE = matrix(rnorm(t*n),t,n);
MU = rep(1,t)%*%t(runif(n,0,4))
X = FF%*%t(W) + EE;

gmm = M2.gmm(X,r = 3,kappa = 0,sigma = rep(1,n),initial = "MLE",eps = 10^-4)
pca = M2.pca(X,r=3,method = "PCA")
ppca = M2.pca(X,C = C,r = 3,method = "P-PCA",J = 4)
TraceRatio(gmm$u,W)
TraceRatio(pca$u,W)
TraceRatio(ppca$u,W)

TraceRatio(gmm$f,FF)
TraceRatio(pca$f,FF)
TraceRatio(ppca$f,FF)

kap_tol = seq(0,5,0.1)
rep = 100
hgl = 1
con_tol = vector()
for (kappa in kap_tol){

  con = vector()
  for (iii in 1:rep) {
    n = 30;t = 100;d = 1;r = 3;
    g1 = function(x){x};
    g2 = function(x){x};
    g3 = function(x){x};
    C = matrix(rnorm(n*d),n,d);W = matrix(NA,n,r);
    W[,1] <- g1(C);W[,2] <- g2(C);W[,3] <- g3(C);
    FF = matrix(rnorm(t*r),t,r);
    EE = matrix(rnorm(t*n),t,n);
    MU = rep(1,t)%*%t(runif(n,2,2))
    X = FF%*%t(W) + EE;

    gmm = M2.gmm(X,r = 3,kappa = kappa,sigma = rep(1,n),initial = "PCA")

    con[iii] <- as.numeric(TraceRatio(gmm$u,W))
  }

  con_tol = cbind(con_tol,con)

  print(hgl)
  hgl = hgl +1
}

plot(colMeans(con_tol),type="b")
boxplot(con_tol)

M3.gmm <- function(X,r,kappa = 0,initial = c("PCA","MLE"),delta = NULL,eps = 10^-6,...){

  n = NCOL(X)
  t = NROW(X)

  m = n^2 + n + 1
  X1 = scale(X,scale = FALSE)
  Mx = cov(X1)

  #PCA
  u_pca = eigen(Mx)$vectors[,1:r]
  f_pca = X1%*%u_pca
  c_pca = f_pca%*%t(u_pca)
  e_pca = X1 - c_pca

  if(initial == "PCA"){
    u_inl_id = u_pca
    sigma_id = colMeans(e_pca^2)
    if(r == 1){
      sk_f = mean(f_pca^3)
    }else{sk_f = colMeans(f_pca^3)}
    sk_id = colMeans(e_pca^3)
  }
  #use MLE to initialize sigma_e
  if(initial == "MLE"){
    B_k = u_pca
    Me_k = diag(diag(cov(e_pca)),n,n)
    theta_k = cbind(B_k,Me_k)
    bias = 1
    iter = 1
    while ( bias > eps ){
      B_k =  theta_k[,1:r]
      Me_k = theta_k[,(r+1):(r+n)]

      Mx_k = B_k%*%t(B_k) + Me_k

      EFF = t(B_k)%*%solve(Mx_k)%*%Mx%*%solve(Mx_k)%*%B_k + diag(1,r,r) - t(B_k)%*%solve(Mx_k)%*%B_k
      EZF = Mx%*%solve(Mx_k)%*%B_k

      B_kk = EZF%*%solve(EFF)
      Me_kk = diag(diag(Mx-B_kk%*%t(B_k)%*%solve(Mx_k)%*%Mx),n,n)

      theta_kk = cbind(B_kk,Me_kk)

      bias = sum((theta_kk-theta_k)^2)

      theta_k = theta_kk
      iter = iter + 1
    }

    u_inl =  theta_k[,1:r]
    Me_inl = theta_k[,(r+1):(r+n)]
    f_inl = t(solve(t(u_inl)%*%diag(1/diag(Me_inl))%*%u_inl)%*%t(u_inl)%*%diag(1/diag(Me_inl))%*%t(X1))
    c_inl = f_inl%*%t(u_inl)
    e_inl = X1 - c_inl

    ev_inl = eigen(t(c_inl)%*%c_inl/t)
    ev_inl_id = ev_inl$values
    u_inl_id = ev_inl$vectors[,1:r]
    f_inl_id = c_inl%*%u_inl_id
    sigma_id = diag(Me_inl)
    if(r == 1){
      sk_f = mean(f_inl_id^3)
    }else{sk_f = colMeans(f_inl_id^3)}
    sk_id = colMeans(e_inl^3)
  }

  sigma_e = sigma_id
  sk_e = sk_id
  #g(x) functions
  #m1
  v1 = colMeans(X)
  #m2
  v2 = t(X)%*%X/t - diag(sigma_e)
  #m3
  v3 = PerformanceAnalytics::M3.MM()

  V = cbind(v1,v2)

  #
  W_sol = matrix(0,m,m)
  for (tt in 1:t) {
    v1_i = X[tt,]
    v2_i = X[tt,]%*%t(X[tt,]) - diag(sigma)
    V_i = cbind(v1_i,v2_i)

    W_sol = W_sol + t(V_i)%*%(diag(n) - u_inl_id%*%t(u_inl_id))%*%V_i
  }
  W_sol = W_sol/t
  eig_w = eigen(W_sol)
  ev_w = as.numeric(eig_w$values[1:(min(n,t))])
  if(is.null(delta)){
    delta = 1/log(t)
  }
  id = which(ev_w > delta)
  evec = matrix(as.numeric(eig_w$vectors[,id]),m,min(n,t))
  WW = (evec[,id])%*%diag(1/ev_w[id])%*%t(evec[,id])

  eig_gmm =  eigen(kappa*t(X)%*%X/t + V%*%WW%*%t(V))

  uv_gmm = as.numeric(matrix(eig_gmm$vectors,n,r))
  u_gmm = matrix(uv_gmm,n,r)
  ev_gmm = (as.numeric(eig_gmm$values))[1:(min(n,t))]

  f_gmm = X%*%u_gmm

  e_gmm = X - f_gmm%*%t(u_gmm)

  con = list(f = f_gmm,u = u_gmm,e = e_gmm, ev = ev_gmm)

  return(con)
}

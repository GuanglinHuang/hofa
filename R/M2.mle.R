#' Estimating the latent factors using maximum likelihood methods
#' @description Estimating the latent factors and factor loadings in high dimensional factor model using maximum likelihood methods based on the covariance or correlation matrix.
#' @param X A matrix or data frame with t rows (samples) and n columns (variables).
#' @param r The number of factors.
#' @param scale logical. If \code{TRUE}, the variance of columns of X are normalized to 1 before factor estimation.
#' @param method Method to use: "\code{ML}", Maximum Likelihood method in Bai and Li (2012); "\code{QML}", Quasi-Maximum Likelihood method in Bai and Li (2013); "\code{ML-GLS}", Bai and Li (2013)'s two step "ML-GLS" method;
#' "\code{ML-ITE}", Bai and Li (2013)'s Iterative "ML-GLS" method; "\code{ML-EM}", Bai and Li (2013)'s EM algorithm for Maximum Likelihood.
#' @param eps The iteration error, default to 10^-8. Available for Maximum Likelihood methods.
#' @param ar.order An integer. Auto regression lag for the idiosyncratic errors in \code{ML-GLS} and \code{ML-ITE}. If it is null, \code{ar.order} will be selected by AIC.
#' @param ... Any other parameters.
#' @return A list of factors, factor loadings and other information, see below.
#' \itemize{
#'   \item{\code{f}}{  Estimated factors.}
#'   \item{\code{u}}{  Estimated factor loadings.}
#'   \item{\code{e}}{  Estimated errors.}
#'   \item{\code{m2e}}{  Diagonal elements of the covariance matrix of errors, only provided in \code{PCA} and \code{ML}.}
#'   \item{\code{rho}}{  Auto regression coefficients of errors, only provided in \code{ML-GLS}, \code{ML-ITE} and \code{ML-EM}.}
#' }
#' @examples
#' n = 100
#' t = 200
#' k = 2
#' par_f = list(rep(1,k),rep(0.8,k),rep(1,k),rep(Inf,k))
#' par_e = list(1,0,2,Inf)
#' rho_f = c(0.5,0.2)
#' par_cove = list(beta = 0.2,J = n/10,rho = 0.2,msig_e = c(1,5))
#' data = hofa.DGP2(n,t,k,par_f,par_e,par_cove,rho_f)$X
#' M2.mle(data,r = 2,method = "ML-EM")


M2.mle = function(X,r,scale = F,method = c("ML","QML","ML-GLS","ML-ITE","ML-EM"),eps = 10^-6,ar.order = 1,...){

  if(method == "ML-EM" && ar.order != 1){
    stop("ML-EM can only be implement under ar.order = 1")
  }
  if(is.null(r)){
    stop("Provide the number of factors")
  }
  n = NCOL(X)
  t = NROW(X)
  X1 = scale(X,scale = scale)
  Mx = cov(X1)*(t-1)/t

  #PCA
  d_pca = eigen(Mx)$values[1:r]/n
  u_pca = eigen(Mx)$vectors[,1:r]%*%diag(sqrt(d_pca),r,r)*sqrt(n)
  f_pca = X1%*%u_pca%*%diag(1/sqrt(d_pca),r,r)/n
  c_pca = f_pca%*%t(u_pca)
  e_pca = X1 - c_pca
  m2e_pca = diag(diag(cov(e_pca)*(t-1)/t),n,n)
  m2_pca = cov(c_pca)*(t-1)/t + m2e_pca

  u_pca = as.matrix(u_pca)
  f_pca = as.matrix(f_pca)
  e_pca = as.matrix(e_pca)

  #ML
  if(method == "ML"||method == "QML"||method == "ML-GLS"||method == "ML-ITE"||method == "ML-EM"){
    B_k = u_pca #initialize
    m2e =  as.matrix(diag(m2e_pca))

    theta_k = hofa::MLE_BL_cpp(r = r,eps = sqrt(eps),Mx = Mx,Bk = B_k,M2E = m2e)

    u_ml =  theta_k[,1:r]
    Me_ml = diag(theta_k[,r+1])
    f_ml = t(solve(t(u_ml)%*%diag(1/diag(Me_ml))%*%u_ml)%*%t(u_ml)%*%diag(1/diag(Me_ml))%*%t(X1))
    e_ml = X1 - f_ml%*%t(u_ml)
    m2_ml = u_ml%*%t(u_ml) + Me_ml

    u_ml = as.matrix(u_ml)
    f_ml = as.matrix(f_ml)
    e_ml = as.matrix(e_ml)

    if(method == "ML"){
      con = list(f = f_ml,u = u_ml,e = e_ml,m2e = diag(Me_ml))
      #print("ML")
    }else{
      #QML
      u_qml = as.matrix(u_ml)
      f_qml = as.matrix(f_ml)
      e_qml = as.matrix(e_ml)
      if(method == "QML"){
        con = list(f = f_qml,u = u_qml,e = e_qml)
        #print("QML")
      }else{
        #ML-GLS
        #step 1 update u
        u_gls = u_qml
        f_gls = f_qml
        e_gls = e_qml
        me_gls = diag(cov(e_gls))
        rho = list()
        for (kk in 1:10){
          for (i in 1:n){
            if(is.null(ar.order)){
              ari = ar(e_gls[,i],order.max = 5)
            }else{ari = ar(e_gls[,i],aic = F,order.max = ar.order)}
            rhoi = ari$ar
            rho[[i]] = rhoi
            pi = ari$order
            f_gls_ar = matrix(NA,t,r)
            for (j in 1:r) {
              f_gls_ar[,j] = lag.mat(f_gls[,j],p = pi)%*%c(1,-rhoi)
            }
            x_gls_ar = lag.mat(X1[,i],p = pi)%*%c(1,-rhoi)
            lmi = lm(x_gls_ar~f_gls_ar-1)
            u_gls[i,] = lmi$coefficients
          }
          #step 2 update f
          f_gls = t(solve(t(u_gls)%*%diag(1/me_gls)%*%u_gls)%*%t(u_gls)%*%diag(1/me_gls)%*%t(X1))
          e_gls = X1 - f_gls%*%t(u_gls)
          me_gls = diag(cov(e_gls))

          u_gls = as.matrix(u_gls)
          f_gls = as.matrix(f_gls)
          e_gls = as.matrix(e_gls)

          if(method == "ML-GLS"||method == "ML-EM"){
            con = list(f = f_gls,u = u_gls,e = e_gls,rho = rho)
            #print("ML-GLS")
            if(method == "ML-EM"){
              f_em = f_gls
              u_em = u_gls
              e_em = e_gls
              rho_em = unlist(rho)
              me_em = diag(cov(e_gls))
              m2u = diag(r)
              error = 1
              iter = 1
              er_tol = vector()
              while (error > eps && iter < 50) {
                Ft = cbind(u_em,-diag(rho_em)%*%u_em)
                Gt = cbind(rbind(diag(r),diag(r)),matrix(0,2*r,r))
                Vt = diag(me_em,n,n)
                Wt = dlm::bdiag(m2u,diag(rep(0,r),r,r))
                Zt = X1[2:t,] - X1[1:(t-1),]%*%diag(rho_em)
                kalman_mod = dlm::dlm(FF = Ft,V = Vt,GG = Gt,W = Wt,m0 = rep(0,2*r), C0 = diag(c(rep(1,r),rep(0,r))))
                kalman_smooth = dlm::dlmSmooth(Zt,mod = kalman_mod)
                f_smooth = as.matrix(kalman_smooth$s)
                V00 = cov(kalman_smooth$s)[1:r,1:r]
                V01 = cov(kalman_smooth$s)[1:r,(r+1):(2*r)]
                V11 = cov(kalman_smooth$s)[(r+1):(2*r),(r+1):(2*r)]
                # update
                u_em_1 = u_em
                rho_em_1 = rho_em
                me_em_1 = me_em
                for (i in 1:n) {
                  u_em_1[i,] = t(solve(V00-rho_em[i]*V01-rho_em[i]*t(V01)+(rho_em[i])^2*V11)%*%as.matrix(apply(as.matrix((f_smooth[2:t,1:r]-rho_em[i]*f_smooth[2:t,(r+1):(2*r)])*(X1[2:t,i] - X1[1:(t-1),i]*rho_em[i])),2,mean)))
                  rho_em_1[i] = sum((X1[2:t,i]*X1[1:(t-1),i]-X1[2:t,i]*as.matrix(f_smooth[2:t,(r+1):(2*r)])%*%u_em_1[i,]-X1[1:(t-1),i]*as.matrix(f_smooth[2:t,(1):(r)])%*%u_em_1[i,]+ rep(u_em_1[i,]%*%V01%*%u_em_1[i,],t-1)))/sum(X1[1:(t-1),i]^2-2*X1[1:(t-1),i]*as.matrix(f_smooth[2:t,(r+1):(2*r)])%*%u_em_1[i,]+ rep(u_em_1[i,]%*%V11%*%u_em_1[i,],t-1))
                  me_em_1[i] = mean((X1[2:t,i] - X1[1:(t-1),i]*rho_em_1[i])^2 - 2*(X1[2:t,i] - X1[1:(t-1),i]*rho_em_1[i])*as.matrix(f_smooth[2:t,(1):(r)])%*%u_em_1[i,] + 2*rho_em_1[i]*(X1[2:t,i] - X1[1:(t-1),i]*rho_em_1[i])*as.matrix(f_smooth[2:t,(r+1):(2*r)])%*%u_em_1[i,] + rep(u_em_1[i,]%*%V00%*%u_em_1[i,],t-1)
                                    - 2*rho_em_1[i]*rep(u_em_1[i,]%*%V01%*%u_em_1[i,],t-1) + (rho_em_1[i])^2*rep(u_em_1[i,]%*%V00%*%u_em_1[i,],t-1))
                  psi_em = V01%*%solve(V11)
                }
                m2u = cov(f_smooth[2:t,1:r] - as.matrix(f_smooth[2:t,(r+1):(2*r)])%*%psi_em)

                error = fnorm(u_em_1-u_em)^2 + fnorm(rho_em_1-rho_em)^2 + fnorm(me_em_1-me_em)^2
                er_tol = c(er_tol,error)
                u_em = u_em_1
                rho_em = rho_em_1
                me_em = me_em_1
                iter = iter+1
              }
              f_em = t(solve(t(u_em)%*%diag(1/me_em)%*%u_em)%*%t(u_em)%*%diag(1/me_em)%*%t(X1))
              e_em = X1 - f_em%*%t(u_em)

              u_em = as.matrix(u_em)
              f_em = as.matrix(f_em)
              e_em = as.matrix(e_em)

              con = list(f = f_em,u = u_em,e = e_em,rho = rho_em)
            }
            break
          }
        }
        if(method == "ML-ITE"){
          con = list(f = f_gls,u = u_gls,e = e_gls,rho = rho)
          #print("ML-ITE")
        }
      }
    }
  }

  return(con)
}


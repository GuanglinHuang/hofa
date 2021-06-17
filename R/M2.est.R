# t = 100
# n = 100
# k = 2
# FF = matrix(rnorm(t*k),t,k)
# W = matrix(rnorm(n*k),n,k)
# E = matrix(rnorm(n*t),t,n)
#
# X = FF%*%t(W)+E
#
# X = read.csv("D:\\OneDrive\\my paper\\GERT\\GERT-JOE\\code\\CODE-LIST\\FRED-MD.csv",row.names = 1)
# X = scale(X)

M2.est = function(X,r,scale = F,method = c("PCA","ML","QML","ML-GLS","ML-ITE","ML-EM","PPCA"),eps = 10^-8,ar.order = 1,...){
if(method == "ML-EM" && ar.order != 1){
  stop("ML-EM can only be implement under ar.order = 1")
}
if(is.null(r)){
    stop("Provide the number of factors")
}
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

if(method == "PCA"){
  con = list(f = f_pca,u = u_pca,e = e_pca,m2e = m2e_pca,m2x = m2_pca)
  #print("PCA")
}

#ML
if(method == "ML"||method == "QML"||method == "ML-GLS"||method == "ML-ITE"||method == "ML-EM"){
  B_k = u_pca #initialize
  Me_k = m2e_pca
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

  u_ml =  theta_k[,1:k]
  Me_ml = theta_k[,(k+1):(k+n)]
  f_ml = t(solve(t(u_ml)%*%diag(1/diag(Me_ml))%*%u_ml)%*%t(u_ml)%*%diag(1/diag(Me_ml))%*%t(X1))
  e_ml = X1 - f_ml%*%t(u_ml)
  m2_ml = u_ml%*%t(u_ml) + Me_ml

  if(method == "ML"){
    con = list(f = f_ml,u = u_ml,e = e_ml,m2e = Me_ml,m2x = m2_ml)
    #print("ML")
  }else{
  #QML
    u_qml = u_ml
    f_qml = f_ml
    e_qml = e_ml
    if(method == "QML"){
      con = list(f = f_qml,u = u_qml,e = e_qml)
      #print("QML")
    }else{
      #ML-GLS
      #step 1 update u
      u_gls = u_qml
      f_gls = f_qml
      e_gls = e_qml
      rho = list()
      for (kk in 1:10){
        for (i in 1:n){
          if(is.null(ar.order)){
            ari = ar(e_gls[,i],order.max = max.order)
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
        f_gls = t(solve(t(u_gls)%*%diag(1/diag(Me_ml))%*%u_gls)%*%t(u_gls)%*%diag(1/diag(Me_ml))%*%t(X1))
        e_gls = X1 - f_gls%*%t(u_gls)
        if(method == "ML-GLS"){
          con = list(f = f_gls,u = u_gls,e = e_gls,rho = rho)
          #print("ML-GLS")
          if(method == "ML-EM"){
          f_em = f_gls
          u_em = u_gls
          e_em = e_gls
          rho_em = unlist(rho)
          #update u
          for (i in 1:n) {

          }
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

if(method == "PPCA"){

}

return(con)
}


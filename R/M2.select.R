#' Determine the number of factors based on the covariance or correlation matrix
#'
#' @param X A matrix or data frame with t rows (samples) and n columns (variables).
#' @param scale logical. If \code{TRUE}, the variance of columns of X are normalized to 1 before factor number test.
#' @param rmax The maximum number of factors.
#' @param method Method to use: "\code{ER}" and "\code{GR}" are Ahn and Horenstein(2013)'s ER and GR estimators;
#' "\code{BN-IC3}","\code{BN-PC3}" and "\code{BIC3}" are Bai and Ng(2002)'s IC3, PC3 and BIC3 estimators, respectively;
#' "\code{ON}" is Onatski(2010)'s estimator; "\code{ACT}" is Fan et al.(2020)'s Adjusted Correlation Thresholding method.
#' @param modified logical. Only available for "\code{ER}" and "\code{GR}". If \code{TRUE}, we use modified estimators which can test zero factor.
#' @param ... Any other parameters.
#' @return The number of factors determined by selected approach and the eigenvalues of the covariance matrix.
#' @examples
#' n = 100
#' t = 200
#' k = 2
#' par_f = list(rep(1,k),rep(0.8,k),rep(1,k),rep(Inf,k))
#' par_e = list(1,0,2,Inf)
#' rho_f = c(0.5,0.2)
#' par_cove = list(beta = 0.2,J = n/10,rho = 0.2,msig_e = c(1,5))
#' data = hofa.DGP2(n,t,k,par_f,par_e,par_cove,rho_f)$X
#' M2.select(data,method = "ER")

M2.select = function(X,scale = F,rmax = 8,
                   method = c("ER","GR","BN-IC3","BN-PC3","BIC3","ON","ACT"),modified = F,...)
{

  #demean
  if(NCOL(X) > NROW(X)){
    X = t(X)
  }
  X = as.matrix(scale(X,scale = scale,center = T))
  n = NCOL(X)
  t = NROW(X)
  kmax = rmax

  m = min(n,t)

  m2 = t(X)%*%X/max(n,t)
  ev2 = (eigen(m2)$values)/(n)

  if(method=="ER"||method=="GR"){

    if(modified == T){
      ev2_0 = sum(ev2)/log(m)
    }else{
      ev2_0 = 0
    }

    ev2_new = c(ev2_0,ev2)
    er2 = ev2_new[1:(kmax+1)]/(ev2_new[2:(kmax+2)])

    V2 = sort(cumsum(c(sort(ev2,decreasing = F),ev2_0)),decreasing = T)

    V2_new = V2

    gr2 = (log(V2_new[1:kmax])-log(V2_new[2:(kmax+1)]))/(log(V2_new[2:(kmax+1)])-log(V2_new[3:(kmax+2)]))

    ER = which(gr2 == max(gr2))-1
    GR = which(er2 == max(er2))-1

    if(method == "ER"){
      FN = ER
    }else{
      FN = GR}
  }

  if(method=="BN-IC3"||method=="BN-PC3"||method=="BIC3"){

    bic3 = ic3 = pc3 = V = p = vector()

    u2 = eigen(m2)$vectors[,1:kmax]

    f_bar = X%*%u2

    sigma_hat = sum((f_bar%*%t(u2) - X)^2)/(n*t)

    for (jjj in 1:kmax) {

      V_j = sum((X%*%u2[,c(1:jjj)]%*%t(u2[,c(1:jjj)]) - X)^2)/(n*t)

      V[jjj] = V_j

      p[jjj] = jjj*sigma_hat*((n+t-jjj)*log(n*t))/(n*t)

      bic3[jjj] = V_j + jjj*sigma_hat*((n+t-jjj)*log(n*t))/(n*t)

      ic3[jjj] = log(V_j) + jjj*(log((n*t)/(n+t))*(n+t)/(n*t))

      pc3[jjj] = V_j + jjj*sigma_hat*(log((n*t)/(n+t))*(n+t)/(n*t))

    }

    BIC3 = which(bic3 == min(bic3))
    IC3 = which(ic3 == min(ic3))
    PC3 = which(pc3 == min(pc3))
    if(method=="BN-IC3"){FN = IC3}
    if(method=="BN-PC3"){FN = PC3}
    if(method=="BIC3"){FN = BIC3}}

  if(method=="ON"){
    j0 =  0
    j  =  kmax + 1
    ed = eigen(m2/max(n,t))$values
    js = 0
    while (abs(j - j0) > 0 ) {
      j0 = j

      ly = ed[j0:(j0+4)]

      lx = c(j0-1,j0,j0+1,j0+2,j0+3)^(2/3)

      lxy = data.frame(ly,lx)

      re = lm(ly~lx,lxy)

      delta = 2*abs(re$coefficients[[2]])

      if(length(which(ed[1:(kmax)]- ed[2:(kmax+1)] > delta)) == 0){
        j = j
        break
      }

      j = max(which(ed[1:(kmax)]- ed[2:(kmax+1)] > delta)) + 1

      js = js + 1
      if(js > 100){
        j = j
        break
      }

    }

    FN = j-1

}

  if(method == "ACT"){
    cor2 = cov2cor(m2)
    lambdaZ = eigen(cor2)$values;
    ev2 = lambdaZ
    DD=NULL; lambdaLY=lambdaZ;
    pp=kmax+2; mz=rep(0,pp); dmz=mz; tmz=mz
    for (kk in 1:pp){
    qu=3/4
    lambdaZ1=lambdaZ[-seq(max(0, 1),kk,1)]; z0=qu*lambdaZ[kk]+(1-qu)*lambdaZ[kk+1]
    ttt=c(lambdaZ1-lambdaZ[kk], z0-lambdaZ[kk]);
    y0=(length(lambdaZ1))/(t-1)
    mz[kk]=-(1-y0)/lambdaZ[kk]+y0*mean(na.omit(1/ttt))
    }

    tempn=(-1/mz)[-1]-1-sqrt(n/(t-1));
    temp1=seq(1,kmax,1);
    temp2=cbind(temp1,tempn[1:kmax])
    ACT = max((temp2[,1][temp2[,2]>0]),0)+1

    FN = min(ACT,kmax)
  }

  return(list(R = FN, eigenvalues = ev2))

}




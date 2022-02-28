#' Determine the number of non-Gaussian and Gaussian factors in hofa
#'
#' @param X A matrix or data frame with t rows (samples) and n columns (variables).
#' @param scale logical. If \code{TRUE}, the variance of columns of X are normalized to 1 before Factor Number test.
#' @param rmax The maximum number of factors.
#' @param method Method to use: "\code{GER4}" and "\code{GGR4}", Lu et al.(2021)'s GER4 and GGR4 estimators; "\code{JJR4}", Jondeau et al.(2018)'s threshold approach.
#' @param modified logical. Only available for "\code{GER4}" and "\code{GGR4}". If \code{TRUE}, we use modified estimators which can test zero factor.
#' @param L An integer. Maximum iteration for JJR4 approach, default to 100.
#' @param ... Any other parameters.
#' @return The number of non-Gaussian and Gaussian factors determined by selected approach.
#' @examples
#' n = 100
#' t = 200
#' k = 2
#' par_f = list(rep(1,k),rep(0.8,k),rep(1,k),rep(Inf,k))
#' par_e = list(1,0,2,Inf)
#' rho_f = c(0.5,0.2)
#' par_cove = list(beta = 0.2,J = n/10,rho = 0.2,msig_e = c(1,5))
#' data = hofa.DGP2(n,t,k,par_f,par_e,par_cove,rho_f)$X
#' M4.select(data,method = "GER4")

M4.select = function(X,scale = F,rmax = 8,method = c("GER4","GGR4","JJR4"),
                     modified = F,L = 100,...)
{

  #demean
  X = as.matrix(scale(X,scale = scale,center = T))
  n = NCOL(X)
  t = NROW(X)
  kmax = rmax

  m = min(n,t)

  if(method == "GER4"||method == "GGR4"){
    m2 = t(X)%*%X/t
    ####cokurtosis
    m4c = C4M(X,m2)
    v4 = eigen(m4c)

    ev4 = v4$values/n^4

    ev4[which(ev4<0)] = 0

    ev4=sqrt(ev4)

    if(modified == T){
      ev4_0 = sum(ev4)/(t/n)
    }else{
      ev4_0 = 0
    }


    ev4_new = c(ev4_0,ev4)
    er4 = ev4_new[1:(kmax+1)]/(ev4_new[2:(kmax+2)])

    if(T %in% is.nan(er4)){
      er4 = er4[-which(is.nan(er4))]
    }

    V4 = sort(cumsum(sort(ev4,decreasing = F)),decreasing = T)

    V4_new = c(sum(ev4_new),V4)

    gr4 = (log(V4_new[1:kmax])-log(V4_new[2:(kmax+1)]))/(log(V4_new[2:(kmax+1)])-log(V4_new[3:(kmax+2)]))

    if(T %in% is.nan(gr4)){
      gr4 = gr4[-which(is.nan(gr4))]
    }

    ger4 = which(er4 == max(er4))-1
    ggr4 = which(gr4 == max(gr4))-1

    # the number of gaussian factors selected by ER and GR
    if(ger4 == 0){u4 = v4$vectors[,ger4]}else{u4 = v4$vectors[,1:ger4]}

    f4 = X%*%u4

    Wg4 = X - f4%*%t(u4)

    vg4 = eigen(cov(Wg4))
    ev2g4 = vg4$values/n

    ev2g4_0 = sum(ev2g4)/log(m)
    ev2g4_new = c(ev2g4_0,ev2g4)
    er2g4 = ev2g4_new[1:(kmax+1)]/(ev2g4_new[2:(kmax+2)])

    if(T %in% is.nan(er2g4)){
      er2g4 = er2g4[-which(is.nan(er2g4))]
    }

    if(ggr4 == 0){u4 = v4$vectors[,ggr4]}else{u4 = v4$vectors[,1:ggr4]}

    f4 = X%*%u4

    Wg4_1 = X - f4%*%t(u4)

    vg4_1 = eigen(cov(Wg4_1))
    ev2g4_1 = vg4_1$values/n

    ev2g4_0_1 = sum(ev2g4_1)/log(m)
    ev2g4_new_1 = c(ev2g4_0_1,ev2g4_1)
    er2g4_1 = ev2g4_new_1[1:(kmax+1)]/(ev2g4_new_1[2:(kmax+2)])

    V2g4_1 = sort(cumsum(sort(ev2g4_1,decreasing = F)),decreasing = T)

    V2g4_new_1 = c(sum(ev2g4_new_1),V2g4_1)

    gr2g4 = (log(V2g4_new_1[1:kmax])-log(V2g4_new_1[2:(kmax+1)]))/(log(V2g4_new_1[2:(kmax+1)])-log(V2g4_new_1[3:(kmax+2)]))

    if(T %in% is.nan(gr2g4)){
      gr2g4 = gr2g4[-which(is.nan(gr2g4))]
    }

    ggr2g4 =  which(gr2g4 == max(gr2g4))-1
    ger2g4 =   which(er2g4 == max(er2g4))-1

    if(method == "GER4"){Rh = ger4}
    if(method == "GGR4"){Rh = ggr4}

    if(method == "GER4"){Rg = ger2g4}
    if(method == "GGR4"){Rg = ggr2g4}

    if(method == "GER4"){R = Rh+Rg}
    if(method == "GGR4"){R = Rh+Rg}

    con = list(R = R,Rh = Rh,Rg = Rg,eigenvalues = ev4)
  }
  if(method == "JJR4"){
    datai = scale(X,scale = T)

    skew_ave = sqrt(mean(as.numeric(PerformanceAnalytics::skewness(datai))^2))
    kurt_ave = sqrt(mean(as.numeric(PerformanceAnalytics::kurtosis(datai))^2))

    start = c(0,2,10)
    sgtpra_tol <- optim(start,function(x){EST_sgt(skew_ave,kurt_ave,x)},method = "L-BFGS-B",lower = c(-0.95,1,3),upper = c(0.95,3,Inf))$par


    simulationcon <- matrix(NA,L,3)

    for (i in 1:L) {
      smu <- matrix(rsgt(n*t,0,1,sgtpra_tol[1],sgtpra_tol[2],sgtpra_tol[3]),t,n)
      smu <- scale(smu)
      lamdajc<- JMCA(smu,kmax)[[1]]

      simulationcon[i,1]<-max(lamdajc)
      simulationcon[i,2]<-median(lamdajc)
      simulationcon[i,3]<-min(lamdajc)

    }

    s_value = colMeans(simulationcon)

    #demean
    Z = scale(X,scale = T)

    JJJC <- JMCA(Z,kmax)

    lamdajc <- JJJC[[1]]

    jjr = length(which(lamdajc > s_value[1]))

    con = list(R = jjr)
  }

  return(con)
}

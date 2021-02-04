#' Determine the number of non-Gaussian and Gaussian factors in HOFA
#'
#' @param X A matrix or data frame with t rows (samples) and n columns (variables).
#' @param scale logical. If \code{TRUE}, the variance of columns of X are normalized to 1 before Factor Number test.
#' @param method Method to use: "\code{ER}", Ahn and Horenstein(2013)'s ER and GR estimators; "\code{GER3}", \code{GER3} and
#' \code{GGR3} estimators; "\code{GER3}", \code{GER4} and \code{GGR4} estimators; "\code{JJR}", Jondeau et al.(2018)'s threshold approach
#' @param rmax The maximum number of factors
#' @param L An integer. Maximum iteration for JJR approach, default to 100
#' @param ... Any other parameters
#' @return The number of non-Gaussian and Gaussian factors determined by selected approach.
#' @examples
#' n = 100;t = 200;k = 2;
#' par_f = c(0.5,0,2,2,Inf,Inf);
#' par_e = c(1,0,2,Inf,0,0,0);
#' rho_ar = c(0.5,0.2);
#' data = hofa.sim(n,t,k,par_f,par_e,rho_ar)$X;
#' GER.sel(data,method = "GER3")

GER.sel = function(X,scale = F,rmax = 8,L = 100,
                   method = c("ER","GER3","GER4","JJR"),...)
  ###input interpretation###
  #X     : a TxN matrix
  #scale : whether the variance of the data need to be normalized
  #kmax  : the maximum number of factors we set
  #method: Ahn and Horenstein(2013)'s ER; GER3; GER4; Jondeau et al.(2018)'s JJR method
  #L     : bootstrap times for JJR's threshold
{

  #demean
  X = as.matrix(scale(X,scale = scale,center = T))
  n = NCOL(X)
  t = NROW(X)
  kmax = rmax

  m = min(n,t)
  mh = min(n,sqrt(t))

  if(method == "ER"){
    m2 = t(X)%*%X/t
    ev2 = (eigen(m2)$values)/(n)
    ev2_0 = sum(ev2)/log(m)
    ev2_new = c(ev2_0,ev2)
    er2 = ev2_new[1:(kmax+1)]/(ev2_new[2:(kmax+2)])

    V2 = sort(cumsum(sort(ev2,decreasing = F)),decreasing = T)

    V2_new = c(sum(ev2_new),V2)

    gr2 = (log(V2_new[1:kmax])-log(V2_new[2:(kmax+1)]))/(log(V2_new[2:(kmax+1)])-log(V2_new[3:(kmax+2)]))

    ggr2 =  which(gr2 == max(gr2))-1
    ger2 =   which(er2 == max(er2))-1

    con = data.frame(ER = ger2,GR = ggr2)

    row.names(con) <- NULL
  }
  if(method == "GER3"){
    m3 = M3M(X)

    v3 = eigen(m3)

    ev3 = v3$values/n^3

    ev3[which(ev3<0)] = 0

    ev3 = sqrt(ev3)

    ev3_0 = sum(ev3)/log(m)
    ev3_new = c(ev3_0,ev3)[1:(kmax+2)]
    er3 = ev3_new[-length(ev3_new)]/(ev3_new[-1])

    if(T %in% is.nan(er3)){
      er3 = er3[-which(is.nan(er3))]
    }

    V3 = sort(cumsum(sort(ev3_new,decreasing = F)),decreasing = T)

    V3_new = V3

    gr3 = (log(V3_new[1:kmax])-log(V3_new[2:(kmax+1)]))/(log(V3_new[2:(kmax+1)])-log(V3_new[3:(kmax+2)]))

    if(T %in% is.nan(gr3)){
      gr3 = gr3[-which(is.nan(gr3))]
    }

    ger3 = which(er3 == max(er3))-1
    ggr3 = which(gr3 == max(gr3))-1

    # second-order factor number selection by M3
    if(ger3 == 0){u3 = v3$vectors[,ger3]}else{u3 = v3$vectors[,1:ger3]}

    f3 = X%*%u3

    Wg3 = X - f3%*%t(u3)

    vg3 = eigen(cov(Wg3))
    ev2g3 = vg3$values/n

    ev2g3_0 = sum(ev2g3)/log(m)
    ev2g3_new = c(ev2g3_0,ev2g3)
    er2g3 = ev2g3_new[1:(kmax+1)]/(ev2g3_new[2:(kmax+2)])

    if(T %in% is.nan(er2g3)){
      er2g3 = er2g3[-which(is.nan(er2g3))]
    }

    if(ggr3 == 0){u3 = v3$vectors[,ggr3]}else{u3 = v3$vectors[,1:ggr3]}

    f3 = X%*%u3

    Wg3_1 = X - f3%*%t(u3)

    vg3_1 = eigen(cov(Wg3_1))
    ev2g3_1 = vg3_1$values/n

    ev2g3_0_1 = sum(ev2g3_1)/log(m)
    ev2g3_new_1 = c(ev2g3_0_1,ev2g3_1)
    er2g3_1 = ev2g3_new_1[1:(kmax+1)]/(ev2g3_new_1[2:(kmax+2)])

    V2g3_1 = sort(cumsum(sort(ev2g3_1,decreasing = F)),decreasing = T)

    V2g3_new_1 = c(sum(ev2g3_new_1),V2g3_1)

    gr2g3 = (log(V2g3_new_1[1:kmax])-log(V2g3_new_1[2:(kmax+1)]))/(log(V2g3_new_1[2:(kmax+1)])-log(V2g3_new_1[3:(kmax+2)]))

    if(T %in% is.nan(gr2g3)){
      gr2g3 = gr2g3[-which(is.nan(gr2g3))]
    }

    ggr2g3 =  which(gr2g3 == max(gr2g3))-1
    ger2g3 =   which(er2g3 == max(er2g3))-1

    con = list(Rh = data.frame(GER3 = ger3,GGR3 = ggr3),Rg = data.frame(ER = ger2g3,GR = ggr2g3))
  }
  if(method == "GER4"){
    m2 = t(X)%*%X/t
    ####cokurtosis
    m4c = C4M(X,m2)
    v4 = eigen(m4c)

    ev4 = v4$values/n^4

    ev4[which(ev4<0)] = 0

    ev4=sqrt(ev4)

    ev4_0 = sum(ev4)/log(m)
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

    con = list(Rh = data.frame(GER4 = ger4,GGR4 = ggr4),Rg = data.frame(ER = ger2g4,GR = ggr2g4))
  }
  if(method == "JJR"){
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

    con = list(Rh = jjr,Rg = 0)
  }
  return(con)
}

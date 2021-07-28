#' Determine the number of non-Gaussian and Gaussian factors in hofa
#'
#' @param X A matrix or data frame with t rows (samples) and n columns (variables).
#' @param scale logical. If \code{TRUE}, the variance of columns of X are normalized to 1 before Factor Number test.
#' @param rmax The maximum number of factors.
#' @param method Method to use: "\code{GER3}" and "\code{GGR3}", Lu et al.(2021)'s GER3 and GGR3 estimators; "\code{JJR3}", Jondeau et al.(2018)'s threshold approach.
#' @param modified logical. Only available for "\code{GER3}" and "\code{GGR3}". If \code{TRUE}, we use modified estimators which can test zero factor.
#' @param L An integer. Maximum iteration for JJR3 approach, default to 100.
#' @param ... Any other parameters.
#' @return The number of non-Gaussian and Gaussian factors determined by selected approach.
#' @examples
#' n = 100
#' t = 200
#' k = 2
#' par_f = c(0.5,0,2,2,Inf,Inf)
#' par_e = c(1,0,2,Inf,0,0,0)
#' rho_ar = c(0.5,0.2)
#' data = hofa.sim(n,t,k,par_f,par_e,rho_ar)$X
#' M3.select(data,method = "GER3")

M3.select = function(X,scale = F,rmax = 8,method = c("GER3","GGR3","JJR3"),
                     modified = F,L = 100,...)
{

  #demean
  X = as.matrix(scale(X,scale = scale,center = T))
  n = NCOL(X)
  t = NROW(X)
  kmax = rmax

  m = min(n,t)

  if(method == "GER3"||method == "GGR3"){
    m3 = M3M(X)/t^2/n^3

    v3 = eigen(m3)

    ev3 = v3$values

    ev3[which(ev3<0)] = 0

    ev3 = sqrt(ev3)

    if(modified == T){
      ev3_0 = sum(ev3)/(t/n)
    }else{
      ev3_0 = 0
    }

    ev3_new = c(ev3_0,ev3)[1:(kmax+2)]
    er3 = ev3_new[-length(ev3_new)]/(ev3_new[-1])

    if(T %in% is.nan(er3)){
      er3 = er3[-which(is.nan(er3))]
    }

    V3 = sort(cumsum(c(sort(ev3,decreasing = F),ev3_0)),decreasing = T)

    V3_new = V3

    gr3 = (log(V3_new[1:kmax])-log(V3_new[2:(kmax+1)]))/(log(V3_new[2:(kmax+1)])-log(V3_new[3:(kmax+2)]))

    if(T %in% is.nan(gr3)){
      gr3 = gr3[-which(is.nan(gr3))]
    }

    ger3 = which(er3 == max(er3))-1
    ggr3 = which(gr3 == max(gr3))-1

    if(method == "GER3"){Rh = ger3}
    if(method == "GGR3"){Rh = ggr3}

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
    ger2g3 =  which(er2g3 == max(er2g3))-1

    if(method == "GER3"){Rg = ger2g3}
    if(method == "GGR3"){Rg = ggr2g3}

    if(method == "GER3"){R = Rh+Rg}
    if(method == "GGR3"){R = Rh+Rg}

    con = list(R = R,Rh = Rh,Rg = Rg,eigenvalues = ev3)
  }
  if(method == "JJR3"){

    datai = scale(X,scale = T)

    skew_ave = sqrt(mean(as.numeric(PerformanceAnalytics::skewness(datai))^2))
    kurt_ave = sqrt(mean(as.numeric(PerformanceAnalytics::kurtosis(datai))^2))

    start = c(0,2,10)
    sgtpra_tol <- optim(start,function(x){EST_sgt(skew_ave,kurt_ave,x)},method = "L-BFGS-B",lower = c(-0.95,1,3),upper = c(0.95,3,Inf))$par


    simulationcon <- matrix(NA,L,3)

    for (i in 1:L) {
      smu <- matrix(rsgt(n*t,0,1,sgtpra_tol[1],sgtpra_tol[2],sgtpra_tol[3]),t,n)
      smu <- scale(smu)
      lamdajc <- JMCA(smu,kmax,gamma = c(0,1,0))[[1]]

      simulationcon[i,1] <- max(lamdajc)
      simulationcon[i,2] <- median(lamdajc)
      simulationcon[i,3] <- min(lamdajc)

    }

    s_value = colMeans(simulationcon)

    #demean
    Z = scale(X,scale = T)

    JJJC <- JMCA(Z,kmax,gamma = c(0,1,0))

    lamdajc <- JJJC[[1]]

    jjr = length(which(lamdajc > s_value[1]))

    con = list(R = jjr)
  }

  return(con)
}




#' Implementation for Alternating Least Square (ALS) algorithm in HOFA
#'
#'
#' @param X A matrix or data frame with t rows (samples) and n columns (variables).
#' @param method Algorithm to use: "\code{HFA3}" to estimate factors by third-order multi-cumulant,
#' "\code{HFA4}" to estimate factors by fourth-order multi-cumulant.
#' @param gamma A weighted vector, a 2x1 vector for \code{HFA3} and a 3x1 vector for \code{HFA4}
#' @param rh The number of non-Gaussian factors
#' @param rg The number of Gaussian factors
#' @param rmax The maximum number of factors
#' @param eps The iteration error, default to 10^-8
#' @param ... Any other parameters
#' @return Estimated factors, factor loadings and commmon components.
#' @examples
#' n = 100;t = 200;k = 2;
#' par_f = c(0.5,0,2,2,Inf,Inf);
#' par_e = c(1,0,2,Inf,0,0,0);
#' rho_ar = c(0.5,0.2);
#' data = HOFA.sim(n,t,k,par_f,par_e,rho_ar)$X;
#' HOFA.als(data,method = "HFA3", gamma = c(1,1),rh = 1, rg = 1)

HOFA.als = function(X,method =c("HFA3","HFA4"),gamma,rh,rg,
               rmax = 8,eps = 10^-8,...){
  n = NCOL(X)
  t = NROW(X)
  kmax = rmax
  if(method == "HFA3"){
    gam0 = gamma[1]
    gam1 = gamma[2]
    if(rh > 0){
      c2_max = X
      m3 = M3M(c2_max)
      m2 = M2M(c2_max)

      eig_j = eigen(gam0*m2/n^2 + gam1*m3/n^3)
      evj_inl = eig_j$values/t^2
      erj = evj_inl[1:kmax]/evj_inl[2:(kmax+1)]
      uh = eig_j$vectors[,1:rh]*sqrt(n)
      fh = c2_max%*%uh/n
      ch = c2_max%*%uh%*%t(uh)/n
      if(rg == 0){
        f = fh; u = uh;c = ch;
      }else{

        eps = eps
        s_g = 1
        s_h = 1
        iter = 0

        while (abs(s_g) > eps || abs(s_h) > eps) {
          Wg = c2_max - fh%*%t(uh)
          m2g = M2M(Wg)/n^2
          eig_m2 = eigen(m2g)
          if(rg == 0){ug = eig_m2$vectors[,rg]*sqrt(n)}else{ug = eig_m2$vectors[,1:rg]*sqrt(n)}

          fg = Wg%*%ug/n

          Wh = c2_max - fg%*%t(ug)
          m2 = M2M(Wh)
          m3 = M3M(Wh)
          mj = gam0*m2/n^2 + gam1*m3/n^3
          eig_mj = eigen(mj)
          evj = eig_mj$values

          if(rh == 0){uh_1 = eig_mj$vectors[,rh]*sqrt(n)}else{uh_1 = eig_mj$vectors[,1:rh]*sqrt(n)}

          fh_1 = c2_max%*%uh_1/n

          Wg = c2_max - fh_1%*%t(uh_1)
          m2g = M2M(Wg)/n^2
          eig_m2 = eigen(m2g)
          if(rg == 0){ug_1 = eig_m2$vectors[,rg]*sqrt(n)}else{ug_1 = eig_m2$vectors[,1:rg]*sqrt(n)}
          fg_1 = Wg%*%ug_1/n

          s_h = fnorm(uh-uh_1)
          s_g = fnorm(ug-ug_1)

          fh = fh_1
          uh = uh_1

          iter = iter + 1

        }

        u = cbind(uh,ug)
        f = c2_max%*%cbind(uh,ug)/n
        c = c2_max%*%u%*%t(u)/n
        colnames(f) <- NULL
      }
    }

    if(rh == 0){
      c2_max = X
      eig_j = eigen(m2/n^2)
      u = eig_j$vectors[,1:rg]*sqrt(n)
      f = c2_max%*%u/n
      c = c2_max%*%u%*%t(u)/n
    }
  }
  if(method == "HFA4"){
    gam0 = gamma[1]
    gam1 = gamma[2]
    gam2 = gamma[3]
    if(rh > 0){
      c2_max = X
      m2 = M2M(c2_max)
      m3 = M3M(c2_max)
      m4 = C4M(c2_max,cov(c2_max))

      eig_j = eigen(gam0*m2/n^2 + gam1*m3/n^3 + gam2*m4/n^4)
      evj_inl = eig_j$values/t^2
      erj = evj_inl[1:kmax]/evj_inl[2:(kmax+1)]
      uh = eig_j$vectors[,1:rh]*sqrt(n)
      fh = c2_max%*%uh/n
      ch = c2_max%*%uh%*%t(uh)/n
      if(rg == 0){
        f = fh;c = ch;u = uh;

      }else{
        eps = eps
        s_g = 1
        s_h = 1
        iter = 0

        while (abs(s_g) > eps || abs(s_h) > eps) {
          Wg = c2_max - fh%*%t(uh)
          m2g = M2M(Wg)/n^2
          eig_m2 = eigen(m2g)
          if(rg == 0){ug = eig_m2$vectors[,rg]*sqrt(n)}else{ug = eig_m2$vectors[,1:rg]*sqrt(n)}

          fg = Wg%*%ug/n

          Wh = c2_max - fg%*%t(ug)
          m2 = M2M(Wh)
          m3 = M3M(Wh)
          m4 = C4M(Wh,t(Wh)%*%(Wh)/t)
          mj = gam0*m2/n^2 + gam1*m3/n^3 + gam2*m4/n^3
          eig_mj = eigen(mj)
          evj = eig_mj$values

          if(rh == 0){uh_1 = eig_mj$vectors[,rh]*sqrt(n)}else{uh_1 = eig_mj$vectors[,1:rh]*sqrt(n)}

          fh_1 = c2_max%*%uh_1/n

          Wg = c2_max - fh_1%*%t(uh_1)
          m2g = M2M(Wg)/n^2
          eig_m2 = eigen(m2g)
          if(rg == 0){ug_1 = eig_m2$vectors[,rg]*sqrt(n)}else{ug_1 = eig_m2$vectors[,1:rg]*sqrt(n)}
          fg_1 = Wg%*%ug_1/n

          s_h = fnorm(uh-uh_1)
          s_g = fnorm(ug-ug_1)

          fh = fh_1
          uh = uh_1

          iter = iter + 1

        }
        u = cbind(uh,ug)
        f = c2_max%*%cbind(uh,ug)/n
        c = c2_max%*%u%*%t(u)/n
        colnames(f) <- NULL
      }


    }

    if(rh == 0){
      c2_max = X
      eig_j = eigen(m2/n^2)
      u = eig_j$vectors[,1:rg]*sqrt(n)
      f = c2_max%*%u/n
      c = c2_max%*%u%*%t(u)/n
    }
  }
  return(list(f = f,u = u,cp = c))
}

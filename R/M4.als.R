#' Fourth-order Alternating Least Square (ALS) algorithm
#' @description Estimating the latent factors using Alternating Least Square (ALS) algorithm based on fourth-order multi-cumulant.
#' @param X A matrix or data frame with t rows (samples) and n columns (variables).
#' @param scale logical. If \code{TRUE}, the variance of columns of X are normalized to 1 before factor estimation.
#' @param gamma A weighted vector, defualt to (0,0,1).
#' @param rh The number of non-Gaussian factors.
#' @param rg The number of Gaussian factors.
#' @param eps The iteration error, default to 10^-8.
#' @param ... Any other parameters.
#' @return Estimated factors, factor loadings and errors.
#' @examples
#' n = 100
#' t = 200
#' k = 2
#' par_f = c(0.5,0,2,2,Inf,Inf)
#' par_e = c(1,0,2,Inf,0,0,0)
#' rho_ar = c(0.5,0.2)
#' data = hofa.sim(n,t,k,par_f,par_e,rho_ar)$X
#' M4.als(data,rh = 1, rg = 1)

M4.als = function(X,scale = FALSE,gamma = NULL,rh,rg,
                                  eps = 10^-8,...){
  n = NCOL(X)
  t = NROW(X)
  X = scale(X,scale = scale)

  if(is.null(gamma)){
    gamma = c(0,0,1)
  }
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
      uh = eig_j$vectors[,1:rh]*sqrt(n)
      fh = c2_max%*%uh/n
      eh = c2_max - fh%*%t(uh)
      if(rg == 0){
        f = fh;u = uh;e = eh;

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
        e = c2_max - f%*%t(u)
        colnames(f) <- NULL
      }


    }

    if(rh == 0){
      c2_max = X
      eig_j = eigen(cov(c2_max)/n)
      evj_inl = eig_j$values
      u = eig_j$vectors[,1:rg]*sqrt(n)
      f = c2_max%*%u/n
      c = c2_max%*%u%*%t(u)/n
      e = c2_max - c
    }

  return(list(f = f,u = u,e = e,ev = evj_inl))
}

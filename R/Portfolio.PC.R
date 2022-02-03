#' Portfolio selection with principal component analysis
#' @description An implementation of the PC portfolio weight optimization through parsimonious higher comoments estimation presented in the paper:
#' Lassance and Vrins (2020) - Portfolio selection with parsimonious higher comoments estimation.
#' @param X A matrix or data frame with t rows (samples) and n columns (variables).
#' @param r An integer, the number of factors. Default to \code{NULL}, several methods exist to select the number of factors.
#' @param rmax An integer, the maximum number of factors. Default to \code{r=10}, only need when \code{r=NULL}.
#' @param fn_sel Factor selection criterion (only need when \code{r=NULL}):"\code{ER}" and "\code{GR}", Ahn and Horenstein(2013)'s ER and GR estimators;
#' "\code{IC3}", Bai and Ng(2002)'s IC3; "\code{ED}", Onatski(2010)'s ED criterion.
#' @param Port_obj The portfolio objective function to be used: Modified Value-at-Risk "\code{MVaR}" and Expected Utility "\code{EU}".
#' @param alpha The confidence level of MVaR (only need when \code{Port_obj="MVaR"}), default to \code{alpha = 0.01}.
#' @param gamma The risk averse parameter of CRRA utility function (only need when \code{Port_obj="EU"}), default to \code{gamma = 10}.
#' @param Adjcov The correction method of covariance matrix: "\code{DNL}", Lodit and Wolf(2018)'s Direct NonLinear shrinkage estimation; "\code{LI}", Lodit and Wolf(2004)'s Linear Identity shrinkage estimation;
#' "\code{NONE}", no correction of the covariance matrix.
#' @param shortselling A logical parameter: whether the portfolio is allowed to shortselling, defalut to \code{TRUE}.
#' @param ... Any other parameters.
#' @return A list contains the optimal portfolio weight, the objective function value, the number of factors, the moments of factors and the moments of epsilon.
#'
#' \itemize{
#'   \item{\code{w}}{  Optimal portfolio weight.}
#'   \item{\code{obj}}{  The series of objective function value, the last one is the optimal objective value.}
#'   \item{\code{r}}{  The number of factors.}
#'   \item{\code{mm_factor}}{The moments of factors.}
#'   \item{\code{mm_eps}}{ The moments of epsilons.}
#'   \item{\code{mm_portfolio}}{ The optimal moments of portfolio.}
#' }
#' @examples
#' data(sp500)
#' Result = Portfolio.PC(sp500,fn_sel = "IC3",Port_obj = "MVaR",Adjcov = "NONE")


Portfolio.PC = function(X,r = NULL,rmax = 10,fn_sel = c("ER","GR","IC3","ED"),
                        Port_obj = c("MVaR","EU"),alpha = 0.01,gamma = 10,Adjcov = c("DNL","LI","NONE"),shortselling = T,...){

  n = NCOL(X)
  t = NROW(X)

  X = scale(X,scale = F)

  if(Adjcov == "DNL"){
    Adj <- function(X){DNLshrink(X)$S_hat}
  }
  if(Adjcov == "LI"){
    Adj <- LIshrink
  }
  if(Adjcov == "NONE"){
    Adj <- cov
  }

  m2 = Adj(X)

  #factor number selection
  if(is.null(r)){
    if(fn_sel == "ER"){
      e2 = eigen(m2)

      ev2 = (e2$values)/n
      er2 = ev2[1:(rmax+1)]/(ev2[2:(rmax+2)])

      r = which(er2 == max(er2))

      if(fn_sel == "GR"){
        V2 = sort(cumsum(sort(ev2,decreasing = F)),decreasing = T)
        gr2 = (log(V2[1:rmax])-log(V2[2:(rmax+1)]))/(log(V2[2:(rmax+1)])-log(V2[3:(rmax+2)]))

        r = which(gr2 == max(gr2))
      }

    }
    if(fn_sel == "IC3"){

      u2 = eigen(m2)$vectors[,1:rmax]
      f_bar = X%*%u2
      sigma_hat = sum((f_bar%*%t(u2) - X)^2)/(n*t)

      ic3 = vector()
      for (jjj in 1:rmax) {
        V_j = sum((X%*%u2[,c(1:jjj)]%*%t(u2[,c(1:jjj)]) - X)^2)/(n*t)
        ic3[jjj] = log(V_j) + jjj*(log((n*t)/(n+t))*(n+t)/(n*t))
      }

      r = which(ic3 == min(ic3))
    }
    if(fn_sel == "ED"){
      j0 =  0
      j  =  rmax + 1
      ed = eigen(m2)$values
      js = 0
      while (abs(j - j0) > 0 ) {
        j0 = j

        ly = ed[j0:(j0+4)]

        lx = c(j0-1,j0,j0+1,j0+2,j0+3)^(2/3)

        lxy = data.frame(ly,lx)

        re = lm(ly~lx,lxy)

        delta = 2*abs(re$coefficients[[2]])

        if(length(which(ed[1:(rmax)]- ed[2:(rmax+1)] > delta)) == 0){
          j = j
          break
        }

        j = max(which(ed[1:(rmax)]- ed[2:(rmax+1)] > delta)) + 1

        js = js + 1

        if(js > 100){
          j = j
          break
        }

      }

      r = j-1
    }
  }


    A_pca = eigen(m2)$vector[,1:r]*sqrt(n) #eigenvalue adjust
    D_pca = eigen(m2)$values[1:r]/n
    f_pca = X%*%A_pca/n

    f = f_pca
    A = A_pca
    eps = X - f%*%t(A)

  m2f = cov(f)
  m3f = PerformanceAnalytics::M3.MM(f)
  c4f = CUM(f)
  mmf = list(m2f,m3f,c4f)

  m2e = diag(Adj(eps))
  m3e = Vm(eps,3)
  m4e = Vm(eps,4)
  mme = list(m2e,m3e,m4e)

  if(shortselling == T){
    lb = -1
  }else{
    lb = 0
  }

  if(Port_obj == "MVaR"){
    fun_obj = function(z){
      mmP = Portfolio.Moments.Mat(z,mm_factor = mmf,mm_eps = mme,A = A);
      obj = Obj.MVaR(mmP,alpha = alpha);
      return(obj)
    }
    rsol = pkgcond::suppress_conditions(Rsolnp::solnp(pars = rep(1/n,n),fun = fun_obj,eqfun = function(x){sum(x)},eqB = 1,LB = rep(lb,n),UB = rep(1,n)))
  }

  if(Port_obj == "EU"){
    fun_obj = function(z){
      mmP = Portfolio.Moments.Mat(z,mm_factor = mmf,mm_eps = mme,A = A)
      obj = Obj.EU(mmP,gamma = gamma);return(obj)
    }
    rsol = pkgcond::suppress_conditions(Rsolnp::solnp(pars = rep(1/n,n),fun = fun_obj,eqfun = function(x){sum(x)},eqB = 1,LB = rep(lb,n),UB = rep(1,n)))
  }
  w.opt = rsol$pars
  w.opt = w.opt/sum(w.opt)
  mmP.opt = Portfolio.Moments.Mat(w.opt,mm_factor = mmf,mm_eps = mme,A = A)

  con = list(w = w.opt,obj = rsol$values, r = r, mm_factor = mmf, mm_eps = mme,mm_portfolio = mmP.opt)

  return(con)
}

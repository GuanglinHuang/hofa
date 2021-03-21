#' Portfolio selection with parsimonious higher comoments estimation
#' Portfolio selection with higher-order comoments estimated by independent components analysis.
#' @param X A matrix or data frame with t rows (samples) and n columns (variables).
#' @param k An integer, the number of factors. Default to \code{NULL}, several methods exist to select the number of factors.
#' @param kmax An integer, the maximum number of factors. Default to \code{k=10}, only need when \code{k=NULL}.
#' @param fn_sel Factor selection criterion (only need when \code{k=NULL}):"\code{ER}" and "\code{GR}", Ahn and Horenstein(2013)'s ER and GR estimators;
#' "\code{IC3}", Bai and Ng(2002)'s IC3; "\code{ED}", Onatski(2010)'s ED criterion.
#' @param Port_obj The portfolio objective function to be used: Modified Value-at-Risk "\code{MVaR}" and Expected Utility "\code{EU}".
#' @param alpha The confidence level of MVaR (only need when \code{Port_obj="MVaR"}), default to \code{alpha = 0.01}.
#' @param gamma The risk averse parameter of CRRA utility function (only need when \code{Port_obj="EU"}), default to \code{gamma = 10}.
#' @param Adjcov The correction method of covariance matrix: "\code{DNL}", Lodit and Wolf(2018)'s Direct NonLinear shrinkage estimation; "\code{LI}", Lodit and Wolf(2004)'s Linear Identity shrinkage estimation.
#' @param shortselling A logical parameter: whether the portfolio is allowed to shortselling, defalut to \code{T}.
#' @return A list contains the optimal portfolio weight \code{w}, the objective function value \code{obj}, the number of factors \code{k}, the moments of factors \code{mm_factor}
#' and the moments of epsilon \code{mm_eps}.
#' @examples
#' data(sp500);
#' Result = Port_ICA(sp500,k = 20,fn_sel = "IC3",Port_obj = "MVaR",Adjcov = "DNL");


Port_ICA = function(X,k = NULL,kmax = 10,fn_sel = c("ER","GR","IC3","ED"),
                     Port_obj = c("MVaR","EU"),alpha = 0.01,gamma = 10,Adjcov = c("DNL","LI"),shortselling = T,...){

  n = NCOL(X)
  t = NROW(X)

  X = scale(X,scale = F)

  if(Adjcov == "DNL"){
    Adj <- function(X){DNLshrink(X)$S_hat}
  }else{
    Adj <- LIshrink
  }

  m2 = Adj(X)

  #factor number selection
  if(is.null(k)){
    if(fn_sel == "ER"){
      e2 = eigen(m2)

      ev2 = (e2$values)/n
      er2 = ev2[1:(kmax+1)]/(ev2[2:(kmax+2)])

      k = which(er2 == max(er2))

      if(fn_sel == "GR"){
        V2 = sort(cumsum(sort(ev2,decreasing = F)),decreasing = T)
        gr2 = (log(V2[1:kmax])-log(V2[2:(kmax+1)]))/(log(V2[2:(kmax+1)])-log(V2[3:(kmax+2)]))

        k = which(gr2 == max(gr2))
      }

    }
    if(fn_sel == "IC3"){

      u2 = eigen(m2)$vectors[,1:kmax]
      f_bar = X%*%u2
      sigma_hat = sum((f_bar%*%t(u2) - X)^2)/(n*t)

      ic3 = vector()
      for (jjj in 1:kmax) {
        V_j = sum((X%*%u2[,c(1:jjj)]%*%t(u2[,c(1:jjj)]) - X)^2)/(n*t)
        ic3[jjj] = log(V_j) + jjj*(log((n*t)/(n+t))*(n+t)/(n*t))
      }

      k = which(ic3 == min(ic3))
    }
    if(fn_sel == "ED"){
      j0 =  0
      j  =  kmax + 1
      ed = eigen(m2)$values
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

      k = j-1
    }
  }

  if(k == 1){
    A = eigen(m2)$vector[,1];
    f = X%*%A;
    eps = X - f%*%t(A);
  }else{
    ICA = ica::icafast(X,nc = k, fun = "kurt")

    f = ICA$S
    A = ICA$M
    eps = X - f%*%t(A)
  }


  m2f = Vm(f,2)
  m3f = Vm(f,3)
  m4f = Vm(f,4)
  mmf = list(m2f,m3f,m4f)

  m2e = diag(Adj(eps))
  m3e = Vm(eps,3)
  m4e = Vm(eps,4)
  mme = list(m2e,m3e,m4e)

  m2_xhat = A%*%t(A) + diag(m2e)

  if(shortselling == T){
    lb = -1
  }else{
    lb = 0
  }

  if(Port_obj == "MVaR"){
    rsol = pkgcond::suppress_conditions(Rsolnp::solnp(pars = rep(1/n,n),fun = function(z){Obj_MVaR(z,mm_factor = mmf,mm_eps = mme,m2_xhat = m2_xhat,A = A,alpha = alpha)},
                         eqfun = function(x){sum(x)},eqB = 1,LB = rep(lb,n),UB = rep(1,n)))
  }

  if(Port_obj == "EU"){
    rsol = pkgcond::suppress_conditions(Rsolnp::solnp(pars = rep(1/n,n),fun = function(z){Obj_EU(z,mm_factor = mmf,mm_eps = mme,m2_xhat = m2_xhat,A = A,gamma = gamma)},
                         eqfun = function(x){sum(x)},eqB = 1,LB = rep(lb,n),UB = rep(1,n)))
  }

  con = list(w = rsol$pars,obj = tail(rsol$values,1), k = k, mm_factor = mmf, mm_eps = mme)

  return(con)
}

Obj_MVaR = function(w,mm_factor,mm_eps,m2_xhat,A,alpha = 0.01,...){
  w = w/sum(w);

  m2f = mm_factor[[1]];
  m3f = mm_factor[[2]];
  m4f = mm_factor[[3]];

  m2e = mm_eps[[1]];
  m3e = mm_eps[[2]];
  m4e = mm_eps[[3]];

  zalpha = qnorm(alpha);
  B = t(w)%*%A
  #M2
  stdP = sqrt(t(w)%*%m2_xhat%*%w);
  #M3
  m3P = sum((B^3)*m3f) + sum((w^3)*m3e);
  skewP = m3P/(stdP^3);
  #M4
  m4P = sum(B^4*(m4f-3*m2f^2)) + sum((w^4)*(m4e-3*m2e^2));
  kurtP = m4P/(stdP^4);

  #MVaR
  Obj = - stdP*zalpha + stdP*(-(1/6)*(zalpha^2-1)*skewP - (1/24)*(zalpha^3-3*zalpha)*kurtP + (1/36)*(2*zalpha^3 - 5*zalpha)*skewP^2);

  return(Obj)
}

Obj_EU = function(w,mm_factor,mm_eps,m2_xhat,A,gamma = 10,...){
  w = w/sum(w);

  m2f = mm_factor[[1]];
  m3f = mm_factor[[2]];
  m4f = mm_factor[[3]];

  m2e = mm_eps[[1]];
  m3e = mm_eps[[2]];
  m4e = mm_eps[[3]];

  B = t(w)%*%A
  #M2
  m2P = t(w)%*%m2_xhat%*%w;
  #M3
  m3P = sum((B^3)*m3f) + sum((w^3)*m3e);
  #M4
  m4P = sum(B^4*(m4f-3*m2f^2)) + sum((w^4)*(m4e-3*m2e^2)) + 3*m2P^2 ;

  #EU
  Obj = - gamma/2*m2P + gamma*(gamma+1)/6*m3P - gamma*(gamma+1)*(gamma+2)/24*m4P

  return(-Obj)
}


Vm = function(X,k){
  X = scale(X,scale = F)
  return(colMeans(X^k))
}

Matmax = function(c = 0,X){
  X[X < c] <- 0
  return(X)
}

DNLshrink = function(X,k = 1,...){
  n = NROW(X)
  p = NCOL(X)
  S = cov(X)*(n-1)/n
  lambda = eigen(S)$values


  ts = which(lambda < 10^(-8))
  if(length(ts)>0){
    lambda[ts] <- lambda[min(ts) - 1]
  }
  u = eigen(S)$vectors


  #compute direct kernel estimator
  lambda = lambda[max(1,p-n+1):p]
  L = lambda%*%t(rep(1,min(p,n)))
  h = n^(-0.35)
  f_tilde = rowMeans(sqrt(Matmax(0,4*t(L)^2*h^2 - (L - t(L))^2))/(2*pi*t(L)^2*h^2))
  Hf_tilde = rowMeans((sign(L - t(L))*sqrt(Matmax(0,(L-t(L))^2 - 4*t(L)^2*h^2))-L+t(L))/
                        (2*pi*t(L)^2*h^2))
  if(p <= n){
    d_tilde = lambda/((pi*(p/n)*lambda*f_tilde)^2 + (1-(p/n)-pi*(p/n)*lambda*Hf_tilde)^2)
  }else{
    Hf_tilde0 = (1-sqrt(1-4*h^2))/(2*pi*h^2)*mean(1/lambda)
    d_tilde0 = 1/(pi*(p-n)/n*Hf_tilde0)
    d_tilde1 = lambda/(pi^2*lambda^2*(f_tilde^2+Hf_tilde^2))
    d_tilde = c(rep(d_tilde0,p-n),d_tilde1)
  }
  d_hat = Iso::pava(d_tilde,decreasing = T)

  S_hat = u%*%diag(d_hat^k)%*%t(u)
  S_sample = u%*%diag(lambda^k)%*%t(u)

  result = list(S_hat = S_hat,
                S_sample = S_sample,
                d_tilde = d_tilde,
                d_hat = d_hat,
                lambda = lambda,
                u      = u,
                f_tilde = f_tilde,
                Hf_tilde = Hf_tilde)
  return(result)
}

LIshrink = function (X, k = 0){
  n <- nrow(X)
  p <- ncol(X)
  if (k == 0) {
    X <- X - tcrossprod(rep(1, n), colMeans(X))
    k = 1
  }
  if (n > k)
    effn <- n - k
  else stop("k must be strictly less than nrow(X)")
  S <- crossprod(X)/effn
  Ip <- diag(p)
  m <- sum(S * Ip)/p
  d2 <- sum((S - m * Ip)^2)/p
  b_bar2 <- 1/(p * effn^2) * sum(apply(X, 1, function(x) sum((tcrossprod(x) -
                                                                S)^2)))
  b2 <- min(d2, b_bar2)
  a2 <- d2 - b2
  return(b2/d2 * m * Ip + a2/d2 * S)
}





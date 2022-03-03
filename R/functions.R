
cov    <- stats::cov
optim  <- stats::optim
rsgt   <- sgt::rsgt
median <- stats::median
rnorm  <- stats::rnorm
runif  <- stats::runif

mm.sgt <- function(h,u,sigma,lambda,p,q)
  ###input interpretation###
  #h: h-th order moment
  #(u,sigma,lambda,p,q): SGT parameters
  {
  mm=vector()
  v=q^(-1/p)*((3*lambda^2+1)*(beta(3/p,q-2/p)/beta(1/p,q))-4*lambda^2*(beta(2/p,q-1/p)/beta(1/p,q))^2)^(-1/2)

  m=(2*v*sigma*lambda*q^(1/p)*beta(2/p,q-1/p))/beta(1/p,q)

  for (r in 0:h) {
    mm[r+1]=choose(h,r)*((1+lambda)^(r+1)+(-1)^(r)*(1-lambda)^(r+1))*(-lambda)^(h-r)*((v*sigma)^h*q^(h/p)/(2^(r-h+1)))*
      (beta((r+1)/p,q-r/p)*beta(2/p,q-1/p)^(h-r)/(beta(1/p,q)^(h-r+1)))
  }
  summm = sum(mm)
  if(is.nan(summm))
    summm = 0
  return(summm)
}

SuperDiag = function(x,k){
  n = length(x)
  if(k==2){
    return(diag(x,n,n))
  }
  if(k == 3){
    mat = matrix(0,n,n^2)
    for (i in 1:n) {
      mat[i,i + (i-1)*n] <- x[i]
    }
    return(mat)
  }
  if(k == 4){
    mat = matrix(0,n,n^3)
    for (i in 1:n) {
      mat[i,i + (i-1)*n + (i-1)*n^2] <- x[i]
    }
    return(mat)
  }
}

CUM <- function(X,...)
  ###input interpretation###
  #X: a TxN matrix
  {
  if(ncol(X) == 0){
    cum <- PerformanceAnalytics::M4.MM(X)
  }else{
    n = length(X[1,])
    t = length(X[,1])
    m2 <- (t-1)/t*cov(X)
    cum<- MC4smaple(X,m2,n,t)}
  return(cum)
}

C4toM4 <- function(C4,m2,...)
{
  n = NCOL(m2)
  M4 <- C4 + t(c(m2))%x%m2 + matrix((m2)%x%(m2),n,n^3) + m2%x%t(c(m2))
  return(M4)
}

M4toC4 <- function(M4,m2,...)
{
  n = NCOL(m2)
  C4 <- M4 - t(c(m2))%x%m2 - matrix((m2)%x%(m2),n,n^3) - m2%x%t(c(m2))
  return(C4)
}

M2M <- function(X)
  ###input interpretation###
  #X: a TxN matrix
  {
  M2M = t(X)%*%((X%*%t(X)))%*%X
  return(M2M)
}

M3M <- function(X)
  ###input interpretation###
  #X: a TxN matrix
  {
  M3M = t(X)%*%((X%*%t(X))*(X%*%t(X)))%*%X
  return(M3M)
}

M4M <- function(X)
  ###input interpretation###
  #X: a TxN matrix
  {
  M4M = t(X)%*%((X%*%t(X))*(X%*%t(X))*(X%*%t(X)))%*%X
  return(M4M)
}

C4M <- function(X,m2)
  ###input interpretation###
  #X: a TxN matrix
  #m2: covariance of X
  {
  t = NROW(X)
  X = scale(X,scale = F)
  M4 = M4M(X)/t^2
  a = apply(X, 1, function(x){sum((x%x%x)*c(m2))})
  b = matrix(a,t,t)
  c = 3*(t(X)%*%((b+t(b))*(X%*%t(X)))%*%X)/t^2
  e1 = 3*as.numeric((t(c(m2))%*%c(m2)))*(m2%*%t(m2))
  e2 = 6*(m2%*%m2%*%m2%*%m2)
  C4 = M4 - c + (e1 + e2)
}


fnorm  <- function(x){sqrt(sum(x^2))}

hofa.lag <- function(X,p){
  c(rep(NA,p),X[1:(length(X)-p)])
}

lag.mat <- function(X,p){
  X.lag = matrix(NA,length(X),p+1)
  for (i in 0:p) {
    X.lag[,i+1] <- hofa.lag(X,i)
  }
  return(X.lag)
}

JMCA = function(X,kmax,gamma = c(1,1,1))
  {
  n  = NROW(X)
  t  = NCOL(X)
  smu  = scale(X,scale = T)
  m2 = M2M(smu)/t^2
  m3 = M3M(smu)/t^2
  m4 = M4M(smu)/t^2

  JJJC <- gamma[1]*m2 + gamma[2]*m3/n^2 + gamma[3]*m4/n^4

  eigval <- sqrt(eigen(JJJC)$values)[1:kmax]

  eigvec <- eigen(JJJC)$vectors[,1:kmax]

  return(list(eigval,eigvec))
}

rccsgt = function(n,ast,p,beta,ut,J,...)
  ###input interpretation###
  #n     : the number of observation T
  #ast   : index for series i
  #p     : the number of response variable N
  #beta,J: control
  #ut    : i.i.d series ut with Tx1
{
  if(ast == 1 )
    et <- ut[,ast]+rowSums(beta*ut[,(ast+1):min(p,(ast+J))])
  if(ast == 2 )
    et <- ut[,ast]+rowSums(beta*ut[,(ast+1):min(p,(ast+J))]) +beta*ut[,1]
  if(ast == (p-1))
    et <- ut[,ast]+rowSums(beta*ut[,(max(1,(ast-J)):(ast-1))]) +beta*ut[,p]
  if(ast == p )
    et <- ut[,ast]+rowSums(beta*ut[,(max(1,(ast-J)):(ast-1))])
  if(ast != p & ast != (p-1) & ast != 1 & ast != 2)
    et <- ut[,ast]+rowSums(beta*ut[,(max(1,(ast-J)):(ast-1))])+rowSums(beta*ut[,(ast+1):min(p,(ast+J))])
  return(et)
}

EST_sgt = function(skew,kurt,theta)
  ###input interpretation###
  #skew  : the average skewness of the dataset
  #kurt  : the average kurtosis of the dataset
  #theta : SGT distribution parameters
{
  lambda =  theta[1]
  p    =  theta[2]
  q    =  theta[3]
  SSE = (mm.sgt(3,0,1,lambda,p,q) - skew)^2 + (mm.sgt(4,0,1,lambda,p,q) - kurt)^2
  return(SSE)
}

TraceRatio = function(f,f0){
  sum(diag(t(f0)%*%f%*%solve(t(f)%*%f)%*%t(f)%*%f0))/sum(diag(t(f0)%*%f0))
}


Portfolio.Moments = function(w,mm_factor,mm_eps,A){

  m2f = mm_factor[[1]];
  m3f = mm_factor[[2]];
  m4f = mm_factor[[3]];

  m2e = mm_eps[[1]];
  m3e = mm_eps[[2]];
  m4e = mm_eps[[3]];

  B = t(w)%*%A

  #M2
  m2P = sum((B^2)*m2f) + sum((w^2)*m2e);
  #M3
  m3P = sum((B^3)*m3f) + sum((w^3)*m3e);
  #M4

  m4P = sum(B^4*(m4f-3*m2f^2)) + sum((w^4)*(m4e-3*m2e^2)) + 3*m2P^2;

  #m4P = sum(B^4*m4f) + 3*sum(t(B^2)%*%(B^2) - diag(diag(t(B^2)%*%(B^2)))) + sum((w^4)*m4e) + 3*sum((w^2*m2e)%*%t(w^2*m2e) - diag(diag((w^2*m2e)%*%t(w^2*m2e)))) + 6*sum((B^2)*m2f)*sum(w^2*m2e)

  mmP = c(m2P,m3P,m4P)
  return(mmP)
}

Portfolio.Moments.GARCHSK = function(w,mm_factor,A){

  m2f = mm_factor[[1]];
  m3f = mm_factor[[2]];
  m4f = mm_factor[[3]];

  B = t(w)%*%A

  #M2
  m2P = sum((B^2)*m2f);
  #M3
  m3P = sum((B^3)*m3f);
  #M4
  m4P = sum((B^4)*m4f);

  mmP = c(m2P,m3P,m4P)
  return(mmP)
}

Portfolio.Moments.Mat = function(w,mm_factor,mm_eps,A){

  m2f = mm_factor[[1]];
  m3f = mm_factor[[2]];
  c4f = mm_factor[[3]];

  m2e = mm_eps[[1]];
  m3e = mm_eps[[2]];
  m4e = mm_eps[[3]];

  B = t(w)%*%A

  #M2
  m2P = B%*%m2f%*%t(B) + sum((w^2)*m2e);
  #M3
  m3P = B%*%m3f%*%(t(B)%x%t(B))+ sum((w^3)*m3e);
  #M4
  m4P = B%*%c4f%*%(t(B)%x%t(B)%x%t(B)) + sum((w^4)*(m4e-3*m2e^2)) + 3*m2P^2;

  mmP = c(m2P,m3P,m4P)
  return(mmP)
}

Obj.MVaR = function(mmP,alpha = 0.01,...){

  m2P = mmP[1];
  m3P = mmP[2];
  m4P = mmP[3];

  zalpha = qnorm(alpha);
  #M2
  stdP = sqrt(m2P);
  #M3
  skewP = m3P/(stdP^3);
  #M4
  kurtP = m4P/(stdP^4);

  #MVaR
  Obj = - stdP*zalpha + stdP*(-(1/6)*(zalpha^2-1)*skewP - (1/24)*(zalpha^3-3*zalpha)*kurtP + (1/36)*(2*zalpha^3 - 5*zalpha)*skewP^2);

  return(Obj)
}

Obj.EU = function(mmP,gamma = 10,...){

  m2P = mmP[1];
  m3P = mmP[2];
  m4P = mmP[3];

  #EU
  Obj =  gamma/2*m2P - gamma*(gamma+1)/6*m3P + gamma*(gamma+1)*(gamma+2)/24*m4P

  return(Obj)
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
  b_bar2 <- 1/(p * effn^2) * sum(apply(X, 1, function(x) sum((tcrossprod(x)-S)^2)))
  b2 <- min(d2, b_bar2)
  a2 <- d2 - b2
  return(b2/d2 * m * Ip + a2/d2 * S)
}

Panel_trans = function(data,type = 1){
  date = row.names(data)
  name = colnames(data)

  n = NCOL(data)
  t = NROW(data)

  panel = vector()
  if(type == 1){
    for (i in 1:n) {
      panel <- rbind(panel,cbind(date,name[i],value[,i]))
    }
    panel_frame = data.frame(date = panel[,1],id = panel[,2],value = as.numeric(panel[,3]))
  }
  if(type == 2){
    for (i in 1:n) {
      panel <- rbind(panel,cbind(name[i],date,value[,i]))
    }
    panel_frame = data.frame(id = panel[,1],date = panel[,2],value = as.numeric(panel[,3]))
  }

  return(panel_frame)
}


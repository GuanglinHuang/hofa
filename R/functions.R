
cov    <- stats::cov
optim  <- stats::optim
rsgt   <- sgt::rsgt
median <- stats::median
rnorm  <- stats::rnorm
runif  <- stats::runif

mm<-function(h,u,sigma,lambda,p,q)
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
    cum<- PerformanceAnalytics::M4.MM(X) - t(c(m2))%x%m2 - matrix((m2)%x%(m2),n,n^3) - m2%x%t(c(m2))}
  return(cum)
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
  M4 = M4M(X)/t^2
  a = apply(X, 1, function(x){sum((x%x%x)*c(m2))})
  b = matrix(a,t,t)
  c = 3*(t(X)%*%((b+t(b))*(X%*%t(X)))%*%X)/t^2
  e1 = 3*as.numeric((t(c(m2))%*%c(m2)))*(m2%*%t(m2))
  e2 = 6*(m2%*%m2%*%m2%*%m2)
  C4 = M4 - c + (e1 + e2)
}

fnorm  <- function(x){sqrt(sum(x^2))}

JMCA = function(X,kmax)
  ###input interpretation###
  #X: a TxN matrix
  #kmax: the maximum number of factors
  {
  n  = NROW(X)
  t  = NCOL(X)
  smu  = scale(X,scale = T)
  m2 = M2M(smu)/t^2
  m3 = M3M(smu)/t^2
  m4 = M4M(smu)/t^2

  JJJC <- m2 + m3/n^2 + m4/n^4

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
  SSE = (mm(3,0,1,lambda,p,q) - skew)^2 + (mm(4,0,1,lambda,p,q) - kurt)^2
  return(SSE)
}



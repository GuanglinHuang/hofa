% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/M2.gmm.R
\name{M2.gmm}
\alias{M2.gmm}
\title{Estimating the latent factors using generalized moment methods}
\usage{
M2.gmm(
  X,
  r,
  kappa = 0,
  sigma_e = NULL,
  initial = c("PCA", "MLE"),
  W_diag = FALSE,
  identity = FALSE,
  delta = NULL,
  eps = 10^-6,
  ...
)
}
\arguments{
\item{X}{A matrix or data frame with t rows (samples) and n columns (variables).}

\item{r}{The number of factors.}

\item{kappa}{An integer. The weight between \code{M} and \code{VWV}, \code{M} denotes the covariance matrix and \code{VWV} denotes the GMM matrix. \code{kappa} default to 0.}

\item{sigma_e}{A \code{n x 1} vector, the variance of the error terms. If it is \code{NULL}, the variance vector will be initialized by PCA or MLE.}

\item{initial}{The method used to initialize the factor loadings and the variance of errors. \code{PCA} denotes Principal Component Analysis and \code{MLE} denotes Maximum Likelihood Estimation in Bai and Li(2012).}

\item{W_diag}{Logical. If \code{TRUE}, the weight matrix \code{W} only has nonzero diagonal elements, default to \code{FALSE}.}

\item{identity}{Logical. If \code{TRUE}, the variances of errors are unified to identity, default to \code{FALSE}.}

\item{delta}{An integer. The hard threshold value of \code{W} matrix in order to guarantee non-singularity of \code{W}, default to 1/log(t).}

\item{eps}{The iteration error, default to 10^-6. Available for initializing the estimators by Maximum Likelihood method.}

\item{...}{Any other parameters.}
}
\value{
A list of factors, factor loadings and other information, see below.
\itemize{
  \item{\code{f}}{  Estimated factors.}
  \item{\code{u}}{  Estimated factor loadings.}
  \item{\code{e}}{  Estimated errors.}
  \item{\code{ev}}{ Eigenvalues of covariance matrix.}
}
}
\description{
Estimating the latent factors and factor loadings in high dimensional factor model using generalized moment methods based on the covariance matrix.
}
\examples{
n = 100
t = 200
k = 2
par_f = list(rep(1,k),rep(0.8,k),rep(1,k),rep(Inf,k))
par_e = list(1,0,2,Inf)
rho_f = c(0.5,0.2)
par_cove = list(beta = 0.2,J = n/10,rho = 0.2,msig_e = c(1,5))
data = hofa.DGP2(n,t,k,par_f,par_e,par_cove,rho_f)$X
M2.gmm(data,r = 2,kappa = 0,sigma_e = rep(1,n),initial = "PCA")
}

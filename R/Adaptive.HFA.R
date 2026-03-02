#' Adaptive Higher-order Factor Analysis (Adaptive.HFA) for weak non-Gaussian factor models
#'
#' @description Estimates latent factors and factor loadings in high-dimensional weak-factor models with non-Gaussianity
#' using an adaptive Higher-order Factor Analysis (HFA) approach. The method adaptively selects optimal cumulant orders
#' (2nd, 3rd, or 4th) based on factor contribution ratios (FCR) to improve estimation accuracy for weak non-Gaussian factors.
#'
#' @param X A matrix or data frame with \code{t} rows (time points/samples) and \code{n} columns (variables).
#' @param scale logical. If \code{TRUE}, columns of \code{X} are scaled to have unit variance before analysis (mean centering is always applied).
#' @param r Integer or NULL. The number of latent factors. If \code{NULL}, the number of factors is automatically determined by
#' combining ER (Eigenvalue Ratio), GER3 (Generalized Eigenvalue Ratio for 3rd cumulant), and GER4 (Generalized Eigenvalue Ratio for 4th cumulant) criteria.
#' @param K Integer (3 or 4 only). Specifies the maximum cumulant order to use for adaptive selection:
#' \itemize{
#'   \item{\code{K=3}: }{Adaptive selection between 2nd (covariance) and 3rd cumulants}
#'   \item{\code{K=4}: }{Adaptive selection between 2nd, 3rd, and 4th cumulants}
#' }
#' @param tau_NT Numeric or NULL. Tuning parameter for factor loading selection. If \code{NULL}, it is automatically set to \code{2*t^(0.25)/n}.
#'
#' @return A list containing the adaptive HFA estimation results, see below:
#' \itemize{
#'   \item{\code{f}}{ Estimated adaptive latent factors (matrix with \code{t} rows and \code{r} columns).}
#'   \item{\code{u}}{ Estimated adaptive factor loadings (matrix with \code{n} rows and \code{r} columns).}
#'   \item{\code{cumulant_order_f}}{ Cumulant order selected for factor estimation (2 = 2nd cumulant/covariance, 3 = 3rd cumulant, 4 = 4th cumulant).}
#'   \item{\code{cumulant_order_u}}{ Cumulant order selected for factor loading estimation (2 = 2nd cumulant/covariance, 3 = 3rd cumulant, 4 = 4th cumulant).}
#' }
#'
#' @details
#' The adaptive HFA method addresses the limitations of traditional PCA (2nd cumulant-based) in weak-factor environments with non-Gaussianity.
#' Key steps include:
#' 1. Preprocessing: Mean centering (always) and optional scaling of the input matrix \code{X}.
#' 2. Factor number determination (if \code{r=NULL}): Combines ER (from covariance) and GER3/GER4 (from higher-order cumulants).
#' 3. Cumulant matrix computation: Calculates 2nd (covariance), 3rd, and 4th cumulant matrices.
#' 4. Factor Contribution Ratio (FCR) calculation: Measures the contribution of top-\code{r} factors to total variation for each cumulant order.
#' 5. Adaptive selection: Chooses the optimal cumulant order for factors/loadings based on maximum FCR (with tuning for loadings).
#'
#' @note
#' - The function relies on auxiliary functions: \code{M2.select}, \code{M3.select}, \code{M4.select} (factor number selection),
#'   \code{M2.pca} (2nd cumulant PCA), \code{M3.als} (3rd cumulant ALS estimation), \code{M4.als} (4th cumulant ALS estimation),
#'   \code{M3M} (3rd cumulant matrix computation), \code{C4M} (4th cumulant matrix computation).
#' - For weak-factor models with pronounced non-Gaussianity (e.g., macroeconomic/financial panels like FRED-MD),
#'   Adaptive.HFA outperforms standard PCA in factor recovery accuracy.
#' @examples
#' n = 100
#' t = 200
#' k = 2
#' par_f = list(rep(1,k),rep(0.8,k),rep(1,k),rep(Inf,k))
#' par_e = list(1,0,2,Inf)
#' rho_f = c(0.5,0.2)
#' par_cove = list(beta = 0.2,J = n/10,rho = 0.2,msig_e = c(1,5))
#' #Non-Gaussian factors
#' data = hofa.DGP2(n,t,k,par_f,par_e,par_cove,rho_f)$X
#' Adaptive.HFA(X = data,scale = FALSE, r = k, K = 3, tau_NT = NULL)
#' #Gaussian factors
#' par_f = list(rep(1,k),rep(0,k),rep(1,k),rep(Inf,k))
#' data = hofa.DGP2(n,t,k,par_f,par_e,par_cove,rho_f)$X
#' Adaptive.HFA(X = data,scale = FALSE, r = k, K = 3, tau_NT = NULL)


Adaptive.HFA = function(X,scale = FALSE, r = NULL, K = 3, tau_NT = NULL){
  X = scale(X,scale = scale)
  n = NCOL(X)
  t = NROW(X)

  if(is.null(tau_NT)){
    tau_NT = 2*t^(0.25)/n
  }

  if(is.null(r)){
    r_er   = M2.select(X,method = "ER", modified = T)$R
    r_ger3 = M3.select(X,method = "GER3")$Rh
    r_ger4 = M4.select(X,method = "GER4")$Rh
    r = max(r_er,r_ger3,r_ger4)
  }

  ev2 = (eigen(cov(X))$values)/(n)
  ev3 = eigen(M3M(X))$values
  ev4 = eigen(C4M(X,cov(X)))$values
  ns = 0.5*n
  FCR2 =  sum(ev2[1:r])/sum(ev2[1:ns])
  FCR3 =  sum(ev3[1:r])/sum(ev3[1:ns])
  FCR4 =  sum(ev4[1:r])/sum(ev4[1:ns])

  res2 = M2.pca(X,r=r,method = "PCA")
  res3 = M3.als(X,gamma = c(0,1),rh = r,rg = 0)
  res4 = M4.als(X,gamma = c(0,0,1),rh = r,rg = 0)

  u2 = res2$u
  u3 = res2$u
  u4 = res2$u

  f2 = res2$f
  f3 = res2$f
  f4 = res2$f

  if (K != 3 & K != 4) {
    stop("Error: The value of K must be either 3 or 4. You provided K = ", K, ".", call. = FALSE)
  }

  if(K == 3){
    place_f = which.max(c(FCR2,FCR3))
    place_u = which.max(c(FCR2+tau_NT,FCR3))

    f_ada = list(f2,f3)[[place_f]]
    u_ada = list(u2,u3)[[place_u]]
  }

  if(K == 4){
    place_f = which.max(c(FCR2,FCR3,FCR4))
    place_u = which.max(c(FCR2+tau_NT,FCR3,FCR4))

    f_ada = list(f2,f3,f4)[[place_f]]
    u_ada = list(u2,u3,u4)[[place_u]]
  }

  return(list(f =f_ada,u = u_ada, cumulant_order_f = place_f + 1, cumulant_order_u = place_u + 1))
}

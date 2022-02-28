#' hofa: Econometric tools for higher-order multi-cumulant factor analysis
#'
#' This R package implements several factor analysis approaches based on the covariance matrix and the higher-order multi-cumulants, including factor number selection, factor estimation and the applications in financial market.
#'
#' @section Factor number selection:
#' The first major part of the project is about determining the number of factors, \code{hofa} implements several approaches based on the covariance matrix and the higher-order moments.
#'
#' The covariance-based approaches: Bai and Ng(2002)'s Information Criterion(IC3,PC3 and BIC3), Onatski(2010)'s Empirical Distribution method(ON), Ahn and Horenstein(2013)'s Eigenvalue Ratio test(ER and GR), Fan et al.(2020)'s Adjusted Correlation Threshold method(ACT). These methods are compiled in \code{M2.select} function.
#'
#' The higher-order moment-based approaches: Lu et al.(2021)'s Generalized Eigenvalue Ratio test(GER3,GER4,GGR3 and GGR4), Jondeau et al.(2018)'s Threshold method(JJR). These methods are compiled in \code{M3.select} and \code{M4.select} functions.
#'
#' @section Factor estimation:
#' The second major part of the project is about factor estimation, \code{hofa} also implements several approaches based on the covariance matrix and the higher-order moments.
#'
#' The covariance-based approaches contain three parts: Principal Component methods, Maximum Likelihood methods and Generalized Moment methods. The \code{M2.pca} function implements classical PCA and Fan et al.(2016)'s Projected PCA(P-PCA). The \code{M2.mle} function implements Bai and Li(2012,2013)'s Maximum Likelihood estimation(ML), Quasi Maximum Likelihood estimation(QML), Generalized Least Square algorithm(ML-GLS), Iterative Generalized Least Square algorithm(ML-ITE) and EM algorithm(ML-EM). The \code{M2.gmm} function implement Fan and Zhong(2018)'s Generalized Moment Method(GMM).
#'
#' The higher-order moment-based approaches: Lu et al.(2021)'s Alternating Least Squares algorithm(\code{M3.als} and \code{M4.als}), Fan and Zhong(2018)'s Generalized Moment Method (\code{M3.gmm}, add third-order moment as structure equations).
#'
#' @section Portfolio selection:
#' The third part of \code{hofa} is about portfolio selection based on the higher-order moments, Lassance and Vrins(2020)'s Independent Component(IC) portfolio and Principal Component(PC) portfolio are implemented in \code{Portfolio.IC} and \code{Portfolio.PC} functions, respectively.
#'
#' @docType package
#' @name hofa
NULL
#> NULL

# [h]igher-[o]rder multi-cumulant [f]actor [a]nalysis

This R package implements several factor analysis approaches based on the covariance matrix and the higher-order multi-cumulants, including factor number selection, factor estimation and the applications in financial market.

Installing the developer version of the package is possible by
```
devtools::install_github("GuanglinHuang/hofa")
```

**Update log:**
Changes in version 0.8.3
 - Remove hofa.als function.
 - Add M3.als function to estimate the factors and factor loadings based on the third-order cumulant.
 - Add M4.als function to estimate the factors and factor loadings based on the fourth-order cumulant.
 
Changes in version 0.8.2
 - Remove M2.est function.
 - Add M2.mle function to estimate the factors based on the covariance by using maximum likelihood methods in Bai and Li(2012,2013).
 - Add M2.pca function to estimate the factors by using principal component analysis (Bai,2003; Fan et al.,2016). 
 
Changes in version 0.8.1
 - Add M2.est function to estimate the factors based on the covariance or correlation matrix. Maximum Likelihood methods in Bai and Li(2012,2013) are available now. This function still works in progress.

Changes in version 0.8.0
 - Remove GER.sel function.
 - Add M2.select function to determine the number of factors based on the second-order moment matrix.
 - Implement Bai and Ng(2002)'s IC3, PC3 and BIC3 estimators in M2.select.
 - Implement Onatski(2010)'s ON estimator in M2.select.
 - Implement Ahn and Horenstein(2013)'s GR and ER estimators in M2.select.
 - Implement Fan et al.(2020)'s ACT estimator in M2.select.
 - Add M3.select and M4.select function to determine the number of factors based on the higher-order moment tensor.
 - Better user parameter setting in M3.select and M4.select functions.

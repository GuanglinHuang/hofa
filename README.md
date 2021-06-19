# [h]igher-[o]rder multi-cumulant [f]actor [a]nalysis

This R package implements several factor analysis approaches based on the covariance matrix and the higher-order multi-cumulants, including factor number selection, factor estimation and the applications in financial market.

Installing the developer version of the package is possible by
```
devtools::install_github("GuanglinHuang/hofa")
```

**Update log:**

Changes in version 0.8.0
 - Remove GER.sel function.
 - Add M2.select function to determine the number of factors based on the second-order moment matrix.
 - Implement Bai and Ng(2002)'s IC3, PC3 and BIC3 estimators in M2.select.
 - Implement Onatski(2010)'s ON estimator in M2.select.
 - Implement Ahn and Horenstein(2013)'s GR and ER estimators in M2.select.
 - Implement Fan et al.(2020)'s ACT estimator in M2.select.
 - Add M3.select and M4.select function to determine the number of factors based on the higher-order moment tensor.
 - Better user parameter setting in M3.select and M4.select functions.
 - a line i wrote on my local computer.

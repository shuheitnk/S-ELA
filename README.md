# Scalarization-based Exploratory Landscape Analysis for Multi-Objective Continuous Optimization Problems

This repository contains code and experiments for Automated Algorithm Selection (AAS) and predicting algorithm performance using our proposed Scalarization-based Exploratory Landscape Analysis (S-ELA).

For now, only the "developer" version is available from GitHub. To install it, either clone this repository or run:

```r
devtools::install_github("shuheitnk/S-ELA")
```

# Quick Start of S-ELA

```r
library(flacco)
library(lhs)
library(smoof)
library(ScalarELA)

d <- 2; n.sample <- 100
fn <- smoof::makeBiObjBBOBFunction(d, fid = 1, iid = 1)

samples <- lhs::improvedLHS(n.sample, d)
X <- sweep(samples, 2, fn.lower <- smoof::getLowerBoxConstraints(fn), "+") * 
     (fn.upper <- smoof::getUpperBoxConstraints(fn) - fn.lower)

Y <- t(apply(X, 1, fn))
Y <- (Y - apply(Y, 2, min)) / (apply(Y, 2, max) - apply(Y, 2, min))
X <- (X - fn.lower) / (fn.upper - fn.lower)

print(DecoELA(X, Y, H = 50, aggregate = TRUE, scalar_func = "weightedsum", set_name = "ela_distr"))
print(DomiELA(X, Y, set_name = "ela_meta"))

```

# Citation

If you are using S-ELA, please use the following BibTeX:

```r
under preparation.
```

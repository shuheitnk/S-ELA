# Scalarization-based Exploratory Landscape Analysis for Multi-Objective Continuous Optimization Problems

This repository contains code and experiments for Automated Algorithm Selection (AAS) and predicting algorithm performance using our proposed Scalarization-based Exploratory Landscape Analysis (S-ELA).

For now, only the "developer" version is available from GitHub. To install it, either clone this repository or run:

```r
devtools::install_github("shuheitnk/S-ELA")
```

# Quick Start of S-ELA

```r
# Load required libraries
library(flacco)
library(lhs)
library(smoof)
library(ScalarELA)

# Set problem dimension and sample size
d <- 2; n.sample <- 100

# Define the bi-objective function from bbob-biobj
fn <- smoof::makeBiObjBBOBFunction(d, fid = 1, iid = 1)

# Generate Latin Hypercube Samples (LHS) in [0,1]^d
samples <- lhs::improvedLHS(n.sample, d)

# Get lower and upper bounds of the function
fn.lower <- smoof::getLowerBoxConstraints(fn)
fn.upper <- smoof::getUpperBoxConstraints(fn)

# Scale samples to match the function's input domain
X <- sweep(samples, 2, fn.lower, "+") * (fn.upper - fn.lower)

# Evaluate the function at the sampled points
Y <- t(apply(X, 1, fn))

# Normalize input variables to [0,1] range
X_scaled <- (X - fn.lower) / (fn.upper - fn.lower)

# Normalize objective values to [0,1] range
Y_min <- apply(Y, 2, min)
Y_max <- apply(Y, 2, max)
Y_scaled <- (Y - Y_min) / (Y_max - Y_min)

# Compute decomposition-based ELA features
print(DecoELA(X_scaled, Y_scaled, H = 50, aggregate = TRUE, scalar_func = "weightedsum", set_name = "ela_distr"))

# Compute non-dominated sorting-based ELA features
print(DomiELA(X_scaled, Y_scaled, set_name = "ela_meta"))


```

# Citation

If you are using S-ELA, please use the following BibTeX:

```r
under preparation.
```

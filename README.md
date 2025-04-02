# Scalarization-based Exploratory Landscape Analysis for Multi-Objective Continuous Optimization Problems

The `ScalarELA` package allows us to implement S-ELA easily.  
For now, only the developer version is available from GitHub. To install it, either clone this repository or run:

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

# Generate sample solutions via Improved Latin Hypercube Sampling in [0,1]^d
samples <- lhs::improvedLHS(n.sample, d)

# Get lower and upper bounds of the decision space
fn.lower <- smoof::getLowerBoxConstraints(fn)
fn.upper <- smoof::getUpperBoxConstraints(fn)

# Scale samples to input the samplued solutions into the function
X <- sweep(samples, 2, fn.lower, "+") * (fn.upper - fn.lower)

# Evaluate the function at the sampled solutions
Y <- t(apply(X, 1, fn))

# Normalize variables to [0,1] range
X_scaled <- (X - fn.lower) / (fn.upper - fn.lower)

# Normalize objective values to [0,1] range
Y_min <- apply(Y, 2, min)
Y_max <- apply(Y, 2, max)
Y_scaled <- (Y - Y_min) / (Y_max - Y_min)

# Compute decomposition-based S-ELA features
print(DecoELA(X_scaled, Y_scaled, H = 5, aggregate = TRUE, scalar_func = "weightedsum", set_name = "ela_distr"))

# Compute NDS-based S-ELA features
print(DomiELA(X_scaled, Y_scaled, set_name = "ela_meta"))


```

# Citation

If you are using S-ELA, please use the following BibTeX:

```r
under preparation.
```

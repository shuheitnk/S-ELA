# ScalarELA

`ScalarELA` is an R package for scalarization-based Exploratory Landscape Analysis (ELA). This package allows you to compute landscape features for various optimization problems.

## Required Packages

To use this package, the following R packages are required:

1. **ScalarELA**: This package itself.
2. **lhs** (version 1.2.0): Provides Latin Hypercube Sampling (LHS) to sample solutions.
3. **smoof** (version 1.6.0.3): Provides bbob-biobj test suite.
4. **MOEADr** (version 1.1.3): Provides Simplex-lattice Design (SLD) for decomposition.
5. **flacco** (version 1.8): Provides the main function to calculate ELA features.

## Installation Instructions

You can install the required packages using the following command:

```r
install.packages(c("lhs", "smoof", "MOEADr", "flacco"))
devtools::install_github("shuheitnk/ScalarELA")
```

Make sure your R version is at least 4.2 (or the compatible version as per package requirements).



# Example Usage

```r
# Load required libraries
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

# Scale samples to input the sampled solutions into the function
X <- sweep(samples, 2, fn.lower, "+") * (fn.upper - fn.lower)

# Evaluate the function at the sampled solutions
Y <- t(apply(X, 1, fn))

# Normalize variables to [0,1] range
X_scaled <- (X - fn.lower) / (fn.upper - fn.lower)

# Normalize objective values to [0,1] range
Y_min <- apply(Y, 2, min)
Y_max <- apply(Y, 2, max)
Y_scaled <- (Y - Y_min) / (Y_max - Y_min)

# Compute decomposition-based S-ELA features (ela_distr)
deco_features = DecoELA(X_scaled, Y_scaled, H = 5, aggregate = TRUE, scalar_func = "weightedsum", set_name = "ela_distr")

# Compute NDS-based S-ELA features (ela_meta)
domi_features = DomiELA(X_scaled, Y_scaled, set_name = "ela_meta")
```

# Citation

If you are using S-ELA, please use the following BibTeX:

```r
under preparation.
```

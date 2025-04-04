# ScalarELA (an R package for implementing S-ELA)
`ScalarELA` is an R package for Scalarization-based Exploratory Landscape Analysis (S-ELA). This package allows you to compute landscape features for various optimization problems. 
S-ELA converts the objective vectors of multi-objective optimization problems into a suitable form, scalar values, for ELA by utilizing decomposition and non-dominated sorting. 
This package has been confirmed to work properly with R version 4.4.1. It is expected to function correctly on this version or later. 
If you are using an earlier version, please ensure R version 3.4.0 or above. 

# Main Functions
## DecoELA
DecoELA convert a multi-objective optimization problem into multiple single-objective sub-problems using weight vectors and apply ELA to each sub-problem. ELA features for each sub-problem are aggregated using descriptive statistics (min, mean, max, sd). Currently, the available scalarization functions are Weighted Sum and Tchebycheff. The decomposition-based approach is implemented by `DecoELA`.
## DomiELA
DomiELA  regard the rank assigned to each solution by non-dominated sorting (NDS) as the objective value and apply ELA. These approaches allow conventional ELA to calculate features of multi-objective continuous optimization problems. The NDS-based approach is implemented by `DomiELA`.

![image](https://github.com/user-attachments/assets/df71a88e-cd1e-44ba-bb1e-da10dba08ffb)


# Installation Instructions
You can install the required packages using the following command:

```r
devtools::install_github("shuheitnk/S-ELA")
```

# Example Usage

```r
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
@inproceedings{,
  title={Scalarization-based Exploratory Landscape Analysis for Multi-Objective Continuous Optimization Problems},
  author={Shuhei, Tanaka and shoichiro, tanaka and toshiharu, hatanaka},
  booktitle={Proceedings of the Genetic and Evolutionary Computation Conference},
  pages={--},
  year={2025}
}
```

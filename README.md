# ScalarELA
`ScalarELA` is an R package for implementing **Scalarization-based Exploratory Landscape Analysis (S-ELA)**. 
This package allows you to compute landscape features for various optimization problems. 

This package has been confirmed to work properly with R version 4.4.1. It is expected to function correctly on this version or later. 
If you are using an earlier version, please ensure R version 3.4.0 or above. 

# Main Functions
## DecoELA
DecoELA convert a multi-objective optimization problem into multiple single-objective sub-problems using weight vectors and apply ELA to each sub-problem. 
ELA features for each sub-problem are aggregated using descriptive statistics (min, mean, max, sd). 
Currently, the available scalarization functions are **Weighted Sum** and **Tchebycheff**. 

The decomposition-based approach is implemented by `DecoELA`.

## DomiELA
DomiELA  regard the rank assigned to each solution by non-dominated sorting (NDS) as the objective value and apply ELA. 
These approaches allow conventional ELA to calculate features of multi-objective continuous optimization problems. 

The NDS-based approach is implemented by `DomiELA`.

![image](https://github.com/user-attachments/assets/df71a88e-cd1e-44ba-bb1e-da10dba08ffb)


# Installation Instructions
You can install `ScalarELA` using the following command:

```r
devtools::install_github("shuheitnk/S-ELA")
```

# Example Usage

```r
library(ScalarELA)

# Set problem dimension and sample size
d <- 2; n.sample <- 100

# Define a bi-objective optimization problem from bbob-biobj
fn <- smoof::makeBiObjBBOBFunction(d, fid = 1, iid = 1)

# Generate sample solutions via Improved Latin Hypercube Sampling in [0,1]^d
samples <- lhs::improvedLHS(n.sample, d)

# Get lower and upper bounds of the decision space
fn.lower <- smoof::getLowerBoxConstraints(fn)
fn.upper <- smoof::getUpperBoxConstraints(fn)

# Scale samples
X <- sweep(samples, 2, fn.lower, "+") * (fn.upper - fn.lower)

# Evaluate the function at sampled solutions
Y <- t(apply(X, 1, fn))

# Compute decomposition-based S-ELA features
deco_features = DecoELA(X,                            # Sampled solutions (matrix)
                        Y,                            # Objective values (matrix)
                        normalize_X = TRUE,           # Whether to normalize sampled solutions
                        normalize_Y = TRUE,           # Whether to normalize objective vectors before decomposition
                        normalize_G = TRUE,           # Whether to normalize objevtive value of sub-problems
                        H = 5,                        # Decomposition parameter H
                        aggregate = TRUE,             # Whether to aggregate sub-problem features
                        scalar_func = "weightedsum",  # Scalarization function (Weighted Sum or Tchebycheff)
                        set_name = "ela_distr"        # Feature set name (ela_meta, ela_distr, disp, nbc, ic, pca, fdc)
                        )

# Compute NDS-based S-ELA features
domi_features = DomiELA(X,
                        Y,
                        normalize_X = TRUE,
                        normalize_R = TRUE,           # Whether to normalize rank values
                        set_name = "ela_meta"
                        )
```

# Citation
This package implements S-ELA which is proposed in the upcoming paper accepted at GECCO2025. 

If you are using S-ELA, please use the following BibTeX:

```r
@inproceedings{,
  title={Scalarization-based Exploratory Landscape Analysis for Multi-Objective Continuous Optimization Problems},
  author={Shuhei, Tanaka and Shoichiro, Tanaka and Toshiharu, Hatanaka},
  booktitle={Proceedings of the Genetic and Evolutionary Computation Conference},
  pages={--},
  year={2025}
}
```

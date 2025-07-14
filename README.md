## ShrinkageTrees
This package provides functions for fitting Horseshoe Trees, Causal Horseshoe Forests, and their more general counterparts, Shrinkage Trees and Causal Shrinkage Forests. These models allow for flexible global-local shrinkage priors on the tree step heights.

The functions can be used for (1) high-dimensional prediction and (2) causal inference of heterogeneous treatment effects given high-dimensional covariates. Outcomes can be continuous, binary, or survival times.

## Installation

You cannot install the released version of ShrinkageTrees from yet:
[CRAN](https://CRAN.R-project.org) with:

``` r
# install.packages("ShrinkageTrees")
```

You can install the development version from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("tijn-jacobs/ShrinkageTrees")
```




## Note to self
 
How to compile the package:
* step 1: (in R) compileAttributes("/Users/tijnjacobs/Library/CloudStorage/OneDrive-VrijeUniversiteitAmsterdam/Documents/GitHub/ShrinkageTrees")
* step 2: (in terminal) R CMD build /Users/tijnjacobs/Library/CloudStorage/OneDrive-VrijeUniversiteitAmsterdam/Documents/GitHub/ShrinkageTrees
* step 3: (in R) install.packages("/Users/tijnjacobs/Library/CloudStorage/OneDrive-VrijeUniversiteitAmsterdam/Documents/GitHub/ShrinkageTrees/ShrinkageTrees_1.0.tar.gz", repos=NULL, type='source')






This package includes data derived from the pdacR package (MIT license), by Richard A. Moffitt. 
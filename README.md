
<!-- README.md is generated from README.Rmd. Please edit that file -->

# MRcML

<!-- badges: start -->

<!-- badges: end -->

Mendelian randomization with constraind maximum likelihood (cML)
methods.

## Installation

<!-- You can install the released version of MRcML from [CRAN](https://CRAN.R-project.org) with:

``` r
install.packages("MRcML")
```
-->

Install the development version from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("xue-hr/MRcML")
```

## Example

Here is an example which shows how to apply MRcML methods to infer the
causal effect from **Fast Glucose (FG)** to **Type-2 Diabetes (T2D)**.

``` r
library(MRcML)
summary(T2D_FG)
#>        Length Class  Mode   
#> b_exp  17     -none- numeric
#> b_out  17     -none- numeric
#> se_exp 17     -none- numeric
#> se_out 17     -none- numeric
```

Example data `T2D_FG` is a list which contains estimated effects sizes
and standard errors of 17 SNPs on T2D and FG. Now we perfrom the main
function with sample size of T2D which is 69033, and using 100 random
start points.

``` r
set.seed(1)
cML_result = mr_cML(T2D_FG$b_exp,
                    T2D_FG$b_out,
                    T2D_FG$se_exp,
                    T2D_FG$se_out,
                    n = 69033,
                    random_start = 100)
```

Now lets take a look at the results:

``` r
cML_result
#> $AIC_theta
#> [1] 1.996619
#> 
#> $AIC_se
#> [1] 0.2194686
#> 
#> $AIC_p
#> [1] 9.242142e-20
#> 
#> $AIC_K
#> [1] 8
#> 
#> $AIC_invalid
#> [1]  1  5  8 12 13 14 15 17
#> 
#> $BIC_theta
#> [1] 2.236579
#> 
#> $BIC_se
#> [1] 0.2061844
#> 
#> $BIC_p
#> [1] 2.050249e-27
#> 
#> $BIC_K
#> [1] 5
#> 
#> $BIC_invalid
#> [1]  8 12 13 15 17
#> 
#> $MA_AIC_theta
#> [1] 2.09175
#> 
#> $MA_AIC_se
#> [1] 0.3012363
#> 
#> $MA_AIC_p
#> [1] 3.814609e-12
#> 
#> $MA_BIC_theta
#> [1] 2.100612
#> 
#> $MA_BIC_se
#> [1] 0.254199
#> 
#> $MA_BIC_p
#> [1] 1.412812e-16
```

BIC selected model gives us indices of invalid IVs: 8, 12, 13, 15, 17.
Now lets draw the scatter plot, invalid IVs are marked with
blue:

<img src="man/figures/README-unnamed-chunk-4-1.png" width="100%" />

<!-- You'll still need to render `README.Rmd` regularly, to keep `README.md` up-to-date.

You can also embed plots, for example:

<img src="man/figures/README-pressure-1.png" width="100%" />

In that case, don't forget to commit and push the resulting figure files, so they display on GitHub! -->

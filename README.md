BGFRA
================
Francisco Javier Luna-Vázquez
2018-02-05

Bayesian Genomic Functional Regression Analysis in R

Instructions for proper implementation
--------------------------------------

### Installation

To complete installation of dev version of BGFRA from GitHub, you have to install a few packages first.

``` r
install.packages('devtools')
devtools::install_github('frahik/BGFRA')
```

### Quick use

``` r
rm(list = ls())
library(BGFRA)
data("wheat_BGFRA")

data <- Wheat # Load from data wheat_BGFRA
Bands <- Bands # Load from data wheat_BGFRA
Wavelengths <- Wavelengths # Load from data wheat_BGFRA

ETA2 <- list(Env = list(X = model.matrix(~0+as.factor(data$Env)), model = "FIXED"),
             Line = list(X = model.matrix(~0+as.factor(data$Line)), model = "BRR"),
             Bands = list(X = Fourier.Basis(Bands, Wavelengths, n.basis = 21, interaction = NULL), model = "BRR")
)

fm2 <- BGFRA(data, ETA = ETA2, nIter = 10000, burnIn= 5000)
```

    ##  Degree of freedom of LP 2  set to default value (5).
    ##  Scale parameter of LP 2  set to default value (4.36140493944161) .
    ##  Degree of freedom of LP 3  set to default value (5).
    ##  Scale parameter of LP 3  set to default value (0.0143843271440224) .

``` r
summary(fm2)
```

    ## --------------------> Summary of data & model <-------------------- 
    ## 
    ##  Number of phenotypes= 300 
    ##  Min (TRN)=  0.3184033 
    ##  Max (TRN)=  7.971726 
    ##  Variance of phenotypes (TRN)= 3.701 
    ##  Residual variance= 0.2529 
    ##  N-TRN= 300   N-TST=0 
    ## 
    ## 
    ##  -- Linear Predictor -- 
    ## 
    ##  Intercept included by default
    ##  Coefficientes in ETA[ 1 ] ( Env ) were asigned a flat prior
    ##  Coefficientes in ETA[ 2 ] ( Line ) modeled as in  BRR 
    ##  Coefficientes in ETA[ 3 ] ( Bands ) modeled as in  BRR 
    ## 
    ## ------------------------------------------------------------------

``` r
pm2 <- BGFRA(data, ETA = ETA2, nIter = 10000, burnIn= 5000, folds = 5, set_seed =10)
```

    ## This might be time demanding, let's take sit and a cup of coffe
    ## 
    ## Done.

``` r
pm2$results
```

    ##           Fold              Env         Cor       MSEP
    ## 1            1        Irrigated  0.20687165 0.14548834
    ## 2            1          Drought  0.22593705 0.44091643
    ## 3            1 ReducedIrrigated  0.61967623 0.06633322
    ## 4            2 ReducedIrrigated  0.39126710 0.18704834
    ## 5            2          Drought  0.62786649 0.14501644
    ## 6            2        Irrigated  0.06380636 0.45363129
    ## 7            3 ReducedIrrigated  0.20638512 0.21356324
    ## 8            3          Drought  0.70661890 0.31135811
    ## 9            3        Irrigated  0.39376149 0.21275294
    ## 10           4 ReducedIrrigated  0.69832098 0.08372290
    ## 11           4        Irrigated -0.29252165 0.47707559
    ## 12           4          Drought  0.83981568 0.17260975
    ## 13           5        Irrigated -0.17877933 0.56578443
    ## 14           5 ReducedIrrigated  0.41622654 0.11252914
    ## 15           5          Drought  0.77003307 0.24211191
    ## 16 Average_all        Irrigated  0.03862771 0.37094652
    ## 17 Average_all          Drought  0.63405424 0.26240253
    ## 18 Average_all ReducedIrrigated  0.46637519 0.13263937

``` r
boxplot(pm2)
```

![](README_files/figure-markdown_github-ascii_identifiers/unnamed-chunk-2-1.png)

### Params

Authors
-------

Francisco Javier Luna-Vázquez (Author, Maintainer)

Osval Antonio Montesinos-López (Author)

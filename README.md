# PAsso

<!-- badges: start -->

[![CRAN checks](https://cranchecks.info/badges/summary/parasol)](https://cran.r-project.org/web/checks/check_results_parasol.html)
[![](http://cranlogs.r-pkg.org/badges/grand-total/parasol?color=blue)](https://cran.r-project.org/package=parasol)
[![](http://cranlogs.r-pkg.org/badges/last-month/parasol?color=green)](https://cran.r-project.org/package=parasol)
[![](http://cranlogs.r-pkg.org/badges/last-week/parasol?color=yellow)](https://cran.r-project.org/package=parasol)
[![](https://travis-ci.org/XiaoruiZhu/parasol.svg?branch=master)](https://travis-ci.org/XiaoruiZhu/parasol)

<!-- badges: end -->

The goal of package PAsso is to assess the Partial Association between ordinal variables.

Overview
--------

An R package of a unified framework for assessing **P**arrtial **Asso**ciation between Ordinal variables. It includes quantification, visualization, and hypothesis testing. All the products are based on the paper by Dungang Liu, Shaobo Li, Yan Yu and Irini Moustaki (2019) and the approach described in [Dungang and Zhang
(2017)](http://www.tandfonline.com/doi/abs/10.1080/01621459.2017.1292915?journalCode=uasa20) and Greenwell et al. (2017, <https://journal.r-project.org/archive/2018/RJ-2018-004/index.html>).

## Installation

The `PAsso` package is currently not available on [PAsso CRAN]() and wait for future updates.

### Install the development version from GitHub

``` r
if (!requireNamespace("devtools")) install.packages("devtools")
devtools::install_github("XiaoruiZhu/PAsso")
```

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(PAsso)
PAsso_1 <- PAsso(responses = c("Prevote.num", "PID"),
                 adjustments = c("income.num", "age", "edu.year"),
                 data = nes2016,
                 uni.model <- "probit"
                 # association = c("partial"),
                 # models = c("probit", "probit"),
                 # method = c("kendall"),
                 # resids.method = "surrogate", fitted.models = NULL,
                 # rep_num = 20
                )

# Print the partial association matrix only
print(PAsso_1, 5)

# Provide partial association matrix, marginal association matrix, and summary of models' coefficients
summary(PAsso_1, 4)

# Plot partial association regression plot: residuals
plot(PAsso_1)

# Association analysis between three ordinal responses ----------
PAsso_2 <- PAsso(responses = c("Prevote.num", "PID", "selfLR"),
                adjustments = c("income.num", "age", "edu.year"),
                data = nes2016,
                uni.model <- "probit",
                method = c("kendall"),
                # models = c("probit", "probit", "probit"),
                # association = c("partial")
                resids.type = "surrogate")

# Compare marginal correlation and partial correlation.
summary(PAsso_2, digits=4)
plot(PAsso_2)

# test function: Conduct inference based on object of "PAsso.test" class ----------------------------
library(progress); #library(doParallel)

system.time(Pcor_SR_test1 <- test(object = PAsso_2, boot_SE=100, H0=0, parallel=F))
print(Pcor_SR_test1, digits=3)

# diagnostic.plot function -----------------------------------------------------
check_qq <- diagnostic.plot(object = PAsso_2, output = "qq")

check_fitted <- diagnostic.plot(object = PAsso_2, output = "fitted")

check_covar <- diagnostic.plot(object = PAsso_2, output = "covariate")

# general association measure and 3-D plot for VOTE and PID ------------------
library("copula")
library("plotly")

testPlots <- plot3D(PAsso_2)
testPlots$plot_1
```

References
----------

Dungang Liu, Shaobo Li, Yan Yu and Irini Moustaki (2019). A unified framework for assessing partial association between ordinal variables: quantification, visualization, and hypothesis testing, working paper.

Dungang Liu & Heping Zhang (2018) Residuals and Diagnostics for Ordinal Regression Models: A Surrogate Approach, Journal of the American Statistical Association, 113:522, 845-854, DOI: 10.1080/01621459.2017.1292915, URL
<http://www.tandfonline.com/doi/abs/10.1080/01621459.2017.1292915?journalCode=uasa20>

Greenwell, B.M., McCarthy, A.J., Boehmke, B.C. & Dungang, L. (2018)
“Residuals and diagnostics for binary and ordinal regression models: An
introduction to the sure package.” The R Journal (pre-print). URL
<https://journal.r-project.org/archive/2018/RJ-2018-004/index.html>


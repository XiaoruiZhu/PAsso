# parasol

<!-- badges: start -->

[![CRAN checks](https://cranchecks.info/badges/summary/parasol)](https://cran.r-project.org/web/checks/check_results_parasol.html)
[![](http://cranlogs.r-pkg.org/badges/grand-total/parasol?color=blue)](https://cran.r-project.org/package=parasol)
[![](http://cranlogs.r-pkg.org/badges/last-month/parasol?color=green)](https://cran.r-project.org/package=parasol)
[![](http://cranlogs.r-pkg.org/badges/last-week/parasol?color=yellow)](https://cran.r-project.org/package=parasol)
[![](https://travis-ci.org/XiaoruiZhu/parasol.svg?branch=master)](https://travis-ci.org/XiaoruiZhu/parasol)

<!-- badges: end -->

The goal of parasol is to assess the Partial Association between Ordinal variables.

Overview
--------

An R package of a unified framework for assessing **Par**rtial **As**sociation between **O**rdina**l** variables. It includes quantification, visualization, and hypothesis testing. All the products are based on the paper by Dungang Liu, Shaobo Li, Yan Yu and Irini Moustaki (2019) and the approach described in [Dungang and Zhang
(2017)](http://www.tandfonline.com/doi/abs/10.1080/01621459.2017.1292915?journalCode=uasa20) and Greenwell et al. (2017, <https://journal.r-project.org/archive/2018/RJ-2018-004/index.html>).

## Installation

The `parasol` package is currently not available on [parasol CRAN]() and wait for future updates.

### Install from GitHub


# Alternatively, install the development version from GitHub

``` r
if (!requireNamespace("devtools")) install.packages("devtools")
devtools::install_github("XiaoruiZhu/parasol")
```

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(parasol)
## Import "nes96" data in "parasol", more details of this data are "?parasol::nes96"
data(nes96)
# Parial association analysis between vote.num and PID:
PAsso_1 <- corr(responses = c("vote.num", "PID"),
                adjustments = c("income.num", "age", "edu.year"),
                data = nes96,
                association = c("partial"),
                models = c("probit", "probit"),
                method = c("kendall"),
                resids.method = "latent",
                fitted.models = NULL, rep_num = 100)

# Marginal association between vote.num and PID:
MAsso_1 <- corr(responses = c("vote.num", "PID"),
                data = nes96,
                association = c("marginal"))

# Compare marginal correlation with partial correlation.
PAsso_1
MAsso_1
```

References
----------

Dungang Liu, Shaobo Li, Yan Yu and Irini Moustaki (2019). A unified framework for assessing partial association between ordinal variables: quantification, visualization, and hypothesis testing.

Liu, D. and Zhang, H. Residuals and Diagnostics for Ordinal Regression
Models: A Surrogate Approach. *Journal of the American Statistical
Association* (accepted). URL
<http://www.tandfonline.com/doi/abs/10.1080/01621459.2017.1292915?journalCode=uasa20>

Greenwell, B.M., McCarthy, A.J., Boehmke, B.C. & Dungang, L. (2018)
“Residuals and diagnostics for binary and ordinal regression models: An
introduction to the sure package.” The R Journal (pre-print). URL
<https://journal.r-project.org/archive/2018/RJ-2018-004/index.html>


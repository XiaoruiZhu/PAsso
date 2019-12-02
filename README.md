# parasol
parasol: An R package for assessing the Partial ASSociation between Ordinal variables
========================================================================================================

Overview
--------

An R package of a unified framework for assessing **Par**rtial **As**sociation between **O**rdina**l** variables. It includes quantification, visualization, and hypothesis testing. All the products are based on the paper [Dungang, Shaobo, Yan (2019)]() and the approach described in [Dungang and Zhang
(2017)](http://www.tandfonline.com/doi/abs/10.1080/01621459.2017.1292915?journalCode=uasa20).

Installation
------------

The `PASSwOrd` package is [currently listed on CRAN]() and can easily be installed:

``` r
# Install from CRAN (recommended)
install.packages("PASSwOrd")

# Alternatively, install the development version from GitHub
if (!requireNamespace("devtools")) install.packages("devtools")
devtools::install_github("XiaoruiZhu/PASSwOrd")
```

References
----------

Dungang, Shaobo, Yan (2019). A nified framework for assessing partial association between ordinal variables: quantification, visualization, and hypothesis testing.

Liu, D. and Zhang, H. Residuals and Diagnostics for Ordinal Regression
Models: A Surrogate Approach. *Journal of the American Statistical
Association* (accepted). URL
<http://www.tandfonline.com/doi/abs/10.1080/01621459.2017.1292915?journalCode=uasa20>

Greenwell, B.M., McCarthy, A.J., Boehmke, B.C. & Dungang, L. (2018)
“Residuals and diagnostics for binary and ordinal regression models: An
introduction to the sure package.” The R Journal (pre-print). URL
<https://journal.r-project.org/archive/2018/RJ-2018-004/index.html>

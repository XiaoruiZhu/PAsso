# PAsso: an R Package for Assessing Partial Association between Ordinal Variables

<!-- badges: start -->

[![](https://img.shields.io/cran/v/PAsso?logo=R)](https://cran.r-project.org/package=PAsso)
[![CRAN checks](https://badges.cranchecks.info/worst/PAsso.svg)](https://cran.r-project.org/web/checks/check_results_PAsso.html)
[![](https://cranlogs.r-pkg.org/badges/grand-total/PAsso?color=blue)](https://cranlogs.r-pkg.org/badges/grand-total/PAsso)
[![](https://cranlogs.r-pkg.org/badges/last-month/PAsso?color=green)](https://cranlogs.r-pkg.org/badges/last-month/PAsso?color=green)
[![](https://cranlogs.r-pkg.org/badges/last-week/PAsso?color=yellow)](https://cranlogs.r-pkg.org/badges/last-week/PAsso?color=yellow)
[![](https://api.travis-ci.com/XiaoruiZhu/PAsso.svg?branch=master)](https://api.travis-ci.com/XiaoruiZhu/PAsso.svg)

<!-- badges: end -->

Overview
--------

An implementation of the unified framework for assessing **P**arrtial **Asso**ciation between ordinal variables after adjusting for a set of covariates (Dungang Liu, Shaobo Li, Yan Yu and Irini Moustaki (2020), accepted by the *Journal of the American Statistical Association*). This package provides a set of tools to quantify, visualize, and test partial associations between multiple ordinal variables. It can produce a number of $\phi$ measures, partial regression plots, 3-D plots, and $p$-values for testing $H_0: \phi=0$ or $H_0: \phi \leq \delta$

## Installation

The `PAsso` package is currently available on [PAsso CRAN](https://CRAN.R-project.org/package=PAsso).

### Install `PAsso` development version from GitHub (recommended)

``` r
# Install the development version from GitHub
if (!requireNamespace("devtools")) install.packages("devtools")
devtools::install_github("XiaoruiZhu/PAsso")
```

### Install `PAsso` from the CRAN

``` r
# Install from CRAN
install.packages("PAsso")

# For macOS, if you have "error: 'math.h' file not found" for installing PAsso v0.1.9,
# the solution could be:
install.packages("https://cran.r-project.org/bin/macosx/el-capitan/contrib/3.6/PAsso_0.1.9.tgz", 
                  repos = NULL, type = "source")
                  
# If error "there is no package called 'gsl'" comes, try:
install.packages("gsl", type = "mac.binary")
```


## Example

The following example shows the R code for evaluating the partial association between a binary variable $\text{PreVote.num}$ and a ordinal variable $\text{PID}$, while adjusting for age, education, and income. Specifically, $\text{PreVote.num}$ is the respondent's voting preference between Donald Trump and Hilary Clinton. And $\text{PID}$ is the respondent’s party identification with 7 ordinal levels from strong democrat (=1) to strong republication (=7). The data set is drawn from the [2016 American National Election Study](https://electionstudies.org/data-center/2016-time-series-study/).

``` r
library(PAsso)

PAsso_1 <- PAsso(responses = c("PreVote.num", "PID"),
                 adjustments = c("income.num", "age", "edu.year"),
                 data = ANES2016,
                 uni.model = "probit",
                 method = c("kendall"))

# Print the partial association matrix only
print(PAsso_1, 5)

# Provide:
# 1. partial association matrix;
# 2. marginal association matrix for comparison purpose;
# 3. summary of models' coefficients for model diagnostics and interpretation
summary(PAsso_1, 4)

# Plot partial association regression plot: residuals
plot(PAsso_1)

# Retrieve residuals that are used as ingredients for partial assocaition analyses
test_resids <- residuals(PAsso_1, draw = 1)
head(test_resids)

# test function: Conduct inference based on object of "PAsso.test" class ----------------------------
library(progress);

system.time(Pcor_SR_test1 <- test(object = PAsso_1, bootstrap_rep = 100, H0 = 0, parallel = F))
print(Pcor_SR_test1, digits=6)

# diagnostic.plot function -----------------------------------------------------
check_qq <- diagnostic.plot(object = PAsso_1, output = "qq")

check_fitted <- diagnostic.plot(object = PAsso_1, output = "fitted")

check_covar <- diagnostic.plot(object = PAsso_1, output = "covariate")

# Or more specific, draw residual-vs-covariate plot for the second model with
# response "PID" and covariate "income.num" 
diagnostic.plot(object = PAsso_1, output = "covariate", x_name = "income.num", model_id = 2)

# general association measure and 3-D plot for VOTE and PID ------------------
library("copula"); library("plotly")

# Draw all pairs
testPlots <- plot3D(PAsso_1)
testPlots$`PreVote.num v.s. PID`

# Draw just one pair
testPlots2 <- plot3D(object = PAsso_1, y1 = "PreVote.num", y2 = "PID")
testPlots2

# "PAsso" advanced using of the function: Input a few models directly ------------------------------
fit_vote <- glm(PreVote.num ~ income.num + age + edu.year, data = ANES2016,
                family = binomial(link = "probit"))
summary(fit_vote)

library(MASS)
fit_PID <- polr(as.factor(PID) ~ income.num + age + edu.year, data = ANES2016,
                method = "probit", Hess = TRUE)

summary(fit_PID)

system.time(PAsso_adv1 <- PAsso(fitted.models=list(fit_vote, fit_PID),
                                association = c("partial"),
                                method = c("kendall"),
                                resids.type = "surrogate")
)

# Partial association coefficients 
print(PAsso_adv1, digits = 3)
summary(PAsso_adv1, digits = 3)
```

References
----------

Li et al. (2021). PAsso: an R Package for Assessing Partial Association between Ordinal Variables. *The R Journal*, 13(2), 239--252, <10.32614/RJ-2021-088>

Liu, D., Li, S., Yu, Y., & Moustaki, I. (2020). Assessing partial association between ordinal variables: quantification, visualization, and hypothesis testing. *Journal of the American Statistical Association*, 1-14. <10.1080/01621459.2020.1796394>

Liu, D., & Zhang, H. (2018). Residuals and diagnostics for ordinal regression models: A surrogate approach. *Journal of the American Statistical Association*, 113(522), 845-854. <10.1080/01621459.2017.1292915>

Greenwell, B.M., McCarthy, A.J., Boehmke, B.C. & Liu, D. (2018) Residuals and diagnostics for binary and ordinal regression models: An
introduction to the sure package. *The R Journal*. <10.32614/RJ-2018-004>

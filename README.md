# PAsso: an R Package for Assessing Partial Association between Ordinal Variables

<!-- badges: start -->

[![](https://www.r-pkg.org/badges/version/PAsso)](https://www.r-pkg.org/badges/version/PAsso)
[![CRAN checks](https://cranchecks.info/badges/summary/PAsso)](https://cran.r-project.org/package=PAsso)
[![](http://cranlogs.r-pkg.org/badges/grand-total/PAsso?color=blue)](http://cranlogs.r-pkg.org/badges/grand-total/PAsso)
[![](http://cranlogs.r-pkg.org/badges/last-month/PAsso?color=green)](http://cranlogs.r-pkg.org/badges/last-month/PAsso?color=green)
[![](http://cranlogs.r-pkg.org/badges/last-week/PAsso?color=yellow)](http://cranlogs.r-pkg.org/badges/last-week/PAsso?color=yellow)
[![](https://travis-ci.com/XiaoruiZhu/PAsso.svg?branch=master)](https://travis-ci.com/XiaoruiZhu/PAsso.svg)

<!-- badges: end -->

Overview
--------

An implementation of the unified framework for assessing **P**arrtial **Asso**ciation between ordinal variables after adjusting for a set of covariates (Dungang Liu, Shaobo Li, Yan Yu and Irini Moustaki (2020), accepted by the *Journal of the American Statistical Association*). This package provides a set of tools to quantify, visualize, and test partial associations between multiple ordinal variables. It can produce a number of $\phi$ measures, partial regression plots, 3-D plots, and $p$-values for testing $H_0: \phi=0$ or $H_0: \phi \leq \delta$

## Installation

The `PAsso` package is currently available on [PAsso CRAN](https://CRAN.R-project.org/package=PAsso).

### Install `PAsso` from the CRAN

``` r
# Install from CRAN
install.packages("PAsso")

# For macOS, if you have "error: 'math.h' file not found" for installing PAsso v0.1.8,
# the solution could be:
install.packages("https://cran.r-project.org/bin/macosx/el-capitan/contrib/3.6/PAsso_0.1.8.tgz", 
                  repos = NULL, type = "source")
                  
# If error "there is no package called 'gsl'" comes, try:
install.packages("gsl", type = "mac.binary")
```

### Install `PAsso` development version from GitHub

``` r
# Alternatively, install the development version from GitHub
if (!requireNamespace("devtools")) install.packages("devtools")
devtools::install_github("XiaoruiZhu/PAsso")
```

## Example

Partial association analysis plays an important role in scientific research. It uncovers the dependency between variables after adjusting for a set of covariates, which are often considered as potential confounding factors. For the data recorded in ordinal scales, the regression models, such as the cumulative link models (McCullagh, 1980), do not have well-defined residuals for assessing the partial association. In the `PAsso` package, we provide a set of tools to quantify, visualize, and test partial associations between multiple ordinal variables.

Specifically, let's consider a pair of ordinal variables $\text{PreVote.num}$ and $\text{PID}$, where the $\text{PreVote.num}$ is the respondent’s voting preference between Donald Trump and Hilary Clinton in the 2016 American National Election Studies. $\text{PID}$ is the respondent’s party identification with 7 ordinal levels from strong democrat (=1) to strong republication (=7). 

One may attempt to answer the following question: **Are the respondent’s voting preference still correlated with their party identification after adjusting for their age, education, and income?** 

To answer this question, we conduct the partial association between $\text{PreVote.num}$ and $\text{PID}$ as an example to illustrate the unified framework implemented in the PAsso package. 

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

Liu, D., Li, S., Yu, Y., & Moustaki, I. (2020) Assessing partial association between ordinal variables: quantification, visualization, and hypothesis testing, accepted by the *Journal of the American Statistical Association*.

Liu, D., & Zhang, H. (2018). Residuals and diagnostics for ordinal regression models: A surrogate approach. *Journal of the American Statistical Association*, 113(522), 845-854. DOI: 10.1080/01621459.2017.1292915, URL
<http://www.tandfonline.com/doi/abs/10.1080/01621459.2017.1292915?journalCode=uasa20>

Greenwell, B.M., McCarthy, A.J., Boehmke, B.C. & Liu, D. (2018)
Residuals and diagnostics for binary and ordinal regression models: An
introduction to the sure package. *The R Journal*. URL <https://journal.r-project.org/archive/2018/RJ-2018-004/index.html>


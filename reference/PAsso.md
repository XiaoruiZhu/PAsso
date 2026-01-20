# Partial association analysis between ordinal responses after adjusting for a set of covariates

This function is mainly designed for conducting the partial association
analysis. It provides two ways of using:

1\. A user-friendly way: only need "responses", "adjustments", and
"data". All the rest of the argument will be setted as default (see
Arguments for details of default).

2\. An advanced way: user can input a list of fitted models by
"fitted.models", then the "responses" and "adjustments" are not
necessary. Supported class of cumulative link models in
[`clm`](https://rdrr.io/pkg/ordinal/man/clm.html),
[`glm`](https://rdrr.io/r/stats/glm.html),
[`lrm`](https://rdrr.io/pkg/rms/man/lrm.html),
[`orm`](https://rdrr.io/pkg/rms/man/orm.html),
[`polr`](https://rdrr.io/pkg/MASS/man/polr.html),
[`vglm`](https://rdrr.io/pkg/VGAM/man/vglm.html), .

It generates an object that has partial association matrix, marginal
association, and some attributes: "arguments" saves c(association,
method, resids.type). "responses" contains the names of response
variables. The attribute "adjustments" contains the names of covariates.
The "summary" function of "PAsso" class of object provides marginal
association ' matrix for comparison purpose.

## Usage

``` r
PAsso(
  responses,
  adjustments,
  data,
  uni.model = c("probit", "logit", "acat"),
  models = NULL,
  method = c("kendall", "pearson", "wolfsigma"),
  resids.type = c("surrogate", "sign", "general", "deviance"),
  jitter = c("latent", "uniform"),
  jitter.uniform.scale = c("probability", "response"),
  fitted.models = NULL,
  n_draws = 20,
  association = "partial",
  ...
)
```

## Arguments

- responses:

  A string vector that specifies response variables. It requires to be
  equal or greater than two variables in the data frame.

- adjustments:

  A string vector specifies covariates/confounders that need to be
  adjusted.

- data:

  A data.frame including responses and adjustments.

- uni.model:

  A character string specifying one single universal model setting for
  all responses. Default `"logit"` refers to cumulative logit model.
  `"probit"` refers to cumulative probit model. `"acat"` fits an
  adjacent categories regression model.

- models:

  A string vector contains default link functions of fitting models with
  respect to each response variable. If `"models"` is missing or has any
  one of the model unspecified, `"uni.model"` is used to specify same
  models for all responses automatically. But, this argument has higher
  priority than the `"uni.model"` as long as the length of `"models"`
  equals to the number of `"responses"`.

- method:

  A string argument to specify correlation coefficient method. Three
  choices `c("kendall", "pearson", "wolfsigma")`. The default is
  `"kendall"`

- resids.type:

  A character string specifying which type of residuals to generate
  Current options include `"surrogate"`, `"sign"`, `"general"`, and
  `"deviance"`. Default is `"surrogate"` residuals.

  `surrogate`

  :   surrogate residuals (Liu and Zhang, 2017);

  `sign`

  :   sign-based residuals (Li and Shepherd, 2010, 2012);

  `general`

  :   generalized residuals (Franses and Paap, 2001);

  `deviance`

  :   deviance residuals (-2\*loglik).

  Although `"sign"`, `"general"`, and `"deviance"` are provided in this
  package, these residuals are problematic for partial association
  analysis between ordinal response (more discussions see Liu, Dungang,
  Li, Shaobo, Yu, Yan, and Moustaki, Irini.(2020))

- jitter:

  A character string specifying how to generate surrogate residuals.
  Current options are `"latent"` and `"uniform"`. Default is `"latent"`.

  `latent`

  :   surrogate residuals.

  `uniform`

  :   sign-based residuals.

- jitter.uniform.scale:

  A character string specifying the scale on which to perform the
  jittering whenever `jitter = "uniform"`. More details:
  [`PAsso::residuals`](http://xiaorui.site/PAsso/reference/residuals.md).

- fitted.models:

  A list contains all the models (S3 objects) you want to assess for the
  partial association between ordinal responses after adjusting for a
  set of covariates covariates. All of these models should be applied to
  the same dataset, having same covariates, same sample size etc. The
  models in this list can be an object of class
  [`clm`](https://rdrr.io/pkg/ordinal/man/clm.html),
  [`glm`](https://rdrr.io/r/stats/glm.html),
  [`lrm`](https://rdrr.io/pkg/rms/man/lrm.html),
  [`orm`](https://rdrr.io/pkg/rms/man/orm.html),
  [`polr`](https://rdrr.io/pkg/MASS/man/polr.html),
  [`vglm`](https://rdrr.io/pkg/VGAM/man/vglm.html).

- n_draws:

  A number to specify draws of surrogate residuls such that the partial
  correlation coefficients are calculated repeatedly. The final
  correlation coefficients are the average of all partial correlation
  coefficients. It is the `"nsim"` argument in `"residuals()"` function.

- association:

  An default argument to specify the partial association. Leave this
  further development of package such that other association analyses
  can be embedded.

- ...:

  Additional optional arguments.

## Value

An object of class `"PAsso"` is a list containing at least the following
components. It contains the partial correlation matrix and multiple
repeats if `n_draws` \> 1. This object has "arguments" attribute saved
as c(association, method, resids.type), "responses" attribute, and
"adjustments" attribute. The list contains:

- `corr`:

  The estimated correlation matrix(average of `rep_MatCorr`) of partial
  association after adjusting confounders;

- `rep_corr`:

  The replications of estimated correlation matrix;

- `rep_SRs`:

  The replications of surrogate residuals if partial association is
  applied;

- `fitted.models`:

  The list stores all fitted.models;

- `data`:

  The data.frame of original dataset;

- `mods_n`:

  The sample size of each fitted model;

- `cor_func`:

  The correlation function after assign different method;

- `Marg_corr`:

  The marginal association matrix.

## References

Liu, D., Li, S., Yu, Y., & Moustaki, I. (2020). Assessing partial
association between ordinal variables: quantification, visualization,
and hypothesis testing. *Journal of the American Statistical
Association*, 1-14.
[doi:10.1080/01621459.2020.1796394](https://doi.org/10.1080/01621459.2020.1796394)

Liu, D., & Zhang, H. (2018). Residuals and diagnostics for ordinal
regression models: A surrogate approach. *Journal of the American
Statistical Association*, 113(522), 845-854.
[doi:10.1080/01621459.2017.1292915](https://doi.org/10.1080/01621459.2017.1292915)

Li, C., & Shepherd, B. E. (2010). Test of association between two
ordinal variables while adjusting for covariates. *Journal of the
American Statistical Association*, 105(490), 612-620.
[doi:10.1198/jasa.2010.tm09386](https://doi.org/10.1198/jasa.2010.tm09386)

Li, C., & Shepherd, B. E. (2012). A new residual for ordinal outcomes.
*Biometrika*, 99(2), 473-480.
[doi:10.1093/biomet/asr073](https://doi.org/10.1093/biomet/asr073)

Franses, P. H., & Paap, R. (2001). *Quantitative models in marketing
research*. Cambridge University Press.
[doi:10.1017/CBO9780511753794](https://doi.org/10.1017/CBO9780511753794)

## Examples

``` r
###########################################################
# User-friendly way of using
###########################################################
library(MASS)

# Import ANES2016 data in "PAsso"
data(ANES2016)

# User-friendly way of the partial association analysis
PAsso_1 <- PAsso(
  responses = c("PreVote.num", "PID"),
  adjustments = c("income.num", "age", "edu.year"),
  data = ANES2016,
  method = c("kendall")
)

print(PAsso_1, digits = 4)
#> -------------------------------------------- 
#> The partial correlation coefficient matrix: 
#>              PreVote.num  PID   
#> PreVote.num  1.0000       0.4499
#> PID                       1.0000
summary(PAsso_1, digits = 4)
#> -------------------------------------------- 
#> The partial correlation coefficient matrix: 
#> 
#>              PreVote.num  PID   
#> PreVote.num  1.0000       0.4499
#> PID                       1.0000
#> -------------------------------------------- 
#> The marginal correlation coefficient matrix: 
#> 
#>              PreVote.num  PID   
#> PreVote.num  1.0000       0.7059
#> PID                       1.0000
#> 
#> --------------------------------------------
#> --------------------------------------------
#> 
#> The coefficients of fitted models are: 
#> 
#>             PreVote.num  PID       
#> income.num   0.0005       0.0009*  
#> Std. Error   0.0005       0.0004   
#> ---                                
#> age          0.0092***    0.0048***
#> Std. Error   0.0016       0.0013   
#> ---                                
#> edu.year    -0.0798***   -0.0459***
#> Std. Error   0.0117       0.0098   
#> ---                                
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

###########################################################
# Advanced way of the partial association analysis
###########################################################

fit.vote <- glm(PreVote.num ~ income.num + age + edu.year,
  data = ANES2016,
  family = binomial(link = "probit")
)
fit.PID <- polr(as.factor(PID) ~ income.num + age + edu.year,
  data = ANES2016,
  method = "probit", Hess = TRUE
)

PAsso_adv1 <- PAsso(
  fitted.models = list(fit.vote, fit.PID),
  method = c("kendall"),
  resids.type = "surrogate"
)

print(PAsso_adv1, digits = 4)
#> -------------------------------------------- 
#> The partial correlation coefficient matrix: 
#>                 PreVote.num  as.factor(PID)
#> PreVote.num     1.0000       0.4493        
#> as.factor(PID)               1.0000        
summary(PAsso_adv1, digits = 4)
#> -------------------------------------------- 
#> The partial correlation coefficient matrix: 
#> 
#>                 PreVote.num  as.factor(PID)
#> PreVote.num     1.0000       0.4493        
#> as.factor(PID)               1.0000        
#> -------------------------------------------- 
#> The marginal correlation coefficient matrix: 
#> 
#>                 PreVote.num  as.factor(PID)
#> PreVote.num     1.0000       0.7059        
#> as.factor(PID)               1.0000        
#> 
#> --------------------------------------------
#> --------------------------------------------
#> 
#> The coefficients of fitted models are: 
#> 
#>             PreVote.num  as.factor(PID)
#> income.num   0.0005       0.0009*      
#> Std. Error   0.0005       0.0004       
#> ---                                    
#> age          0.0092***    0.0048***    
#> Std. Error   0.0016       0.0013       
#> ---                                    
#> edu.year    -0.0798***   -0.0459***    
#> Std. Error   0.0117       0.0098       
#> ---                                    
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
```

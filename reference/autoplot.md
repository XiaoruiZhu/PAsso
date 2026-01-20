# Residual-based diagnostic plots

Residual-based diagnostic plots for cumulative link and general
regression models using
[`ggplot`](https://ggplot2.tidyverse.org/reference/ggplot.html)
graphics.

## Usage

``` r
# S3 method for class 'resid'
autoplot(
  object,
  output = c("qq", "fitted", "covariate"),
  x = NULL,
  fit = NULL,
  distribution = qnorm,
  ncol = NULL,
  alpha = 1,
  xlab = NULL,
  color = "#444444",
  shape = 19,
  size = 2,
  qqpoint.color = "#444444",
  qqpoint.shape = 19,
  qqpoint.size = 2,
  qqline.color = "#888888",
  qqline.linetype = "dashed",
  qqline.size = 1,
  smooth = TRUE,
  smooth.color = "red",
  smooth.linetype = 1,
  smooth.size = 1,
  fill = NULL,
  resp_name = NULL,
  ...
)

# S3 method for class 'glm'
autoplot(
  object,
  output = c("qq", "fitted", "covariate"),
  x = NULL,
  fit = NULL,
  distribution = qnorm,
  ncol = NULL,
  alpha = 1,
  xlab = NULL,
  color = "#444444",
  shape = 19,
  size = 2,
  qqpoint.color = "#444444",
  qqpoint.shape = 19,
  qqpoint.size = 2,
  qqline.color = "#888888",
  qqline.linetype = "dashed",
  qqline.size = 1,
  smooth = TRUE,
  smooth.color = "red",
  smooth.linetype = 1,
  smooth.size = 1,
  fill = NULL,
  resp_name = NULL,
  ...
)

# S3 method for class 'clm'
autoplot(
  object,
  output = c("qq", "fitted", "covariate"),
  x = NULL,
  fit = NULL,
  distribution = qnorm,
  ncol = NULL,
  alpha = 1,
  xlab = NULL,
  color = "#444444",
  shape = 19,
  size = 2,
  qqpoint.color = "#444444",
  qqpoint.shape = 19,
  qqpoint.size = 2,
  qqline.color = "#888888",
  qqline.linetype = "dashed",
  qqline.size = 1,
  smooth = TRUE,
  smooth.color = "red",
  smooth.linetype = 1,
  smooth.size = 1,
  fill = NULL,
  resp_name = NULL,
  ...
)

# S3 method for class 'lrm'
autoplot(
  object,
  output = c("qq", "fitted", "covariate"),
  x = NULL,
  fit = NULL,
  distribution = qnorm,
  ncol = NULL,
  alpha = 1,
  xlab = NULL,
  color = "#444444",
  shape = 19,
  size = 2,
  qqpoint.color = "#444444",
  qqpoint.shape = 19,
  qqpoint.size = 2,
  qqline.color = "#888888",
  qqline.linetype = "dashed",
  qqline.size = 1,
  smooth = TRUE,
  smooth.color = "red",
  smooth.linetype = 1,
  smooth.size = 1,
  fill = NULL,
  resp_name = NULL,
  ...
)

# S3 method for class 'orm'
autoplot(
  object,
  output = c("qq", "fitted", "covariate"),
  x = NULL,
  fit = NULL,
  distribution = qnorm,
  ncol = NULL,
  alpha = 1,
  xlab = NULL,
  color = "#444444",
  shape = 19,
  size = 2,
  qqpoint.color = "#444444",
  qqpoint.shape = 19,
  qqpoint.size = 2,
  qqline.color = "#888888",
  qqline.linetype = "dashed",
  qqline.size = 1,
  smooth = TRUE,
  smooth.color = "red",
  smooth.linetype = 1,
  smooth.size = 1,
  fill = NULL,
  resp_name = NULL,
  ...
)

# S3 method for class 'polr'
autoplot(
  object,
  output = c("qq", "fitted", "covariate"),
  x = NULL,
  fit = NULL,
  distribution = qnorm,
  ncol = NULL,
  alpha = 1,
  xlab = NULL,
  color = "#444444",
  shape = 19,
  size = 2,
  qqpoint.color = "#444444",
  qqpoint.shape = 19,
  qqpoint.size = 2,
  qqline.color = "#888888",
  qqline.linetype = "dashed",
  qqline.size = 1,
  smooth = TRUE,
  smooth.color = "red",
  smooth.linetype = 1,
  smooth.size = 1,
  fill = NULL,
  resp_name = NULL,
  ...
)

# S3 method for class 'vglm'
autoplot(
  object,
  output = c("qq", "fitted", "covariate"),
  x = NULL,
  fit = NULL,
  distribution = qnorm,
  ncol = NULL,
  alpha = 1,
  xlab = NULL,
  color = "#444444",
  shape = 19,
  size = 2,
  qqpoint.color = "#444444",
  qqpoint.shape = 19,
  qqpoint.size = 2,
  qqline.color = "#888888",
  qqline.linetype = "dashed",
  qqline.size = 1,
  smooth = TRUE,
  smooth.color = "red",
  smooth.linetype = 1,
  smooth.size = 1,
  fill = NULL,
  resp_name = NULL,
  ...
)
```

## Arguments

- object:

  An object of class `sure:resids`,
  [`clm`](https://rdrr.io/pkg/ordinal/man/clm.html),
  [`glm`](https://rdrr.io/r/stats/glm.html),
  [`lrm`](https://rdrr.io/pkg/rms/man/lrm.html),
  [`orm`](https://rdrr.io/pkg/rms/man/orm.html),
  [`polr`](https://rdrr.io/pkg/MASS/man/polr.html), or
  [`vglm`](https://rdrr.io/pkg/VGAM/man/vglm.html).

- output:

  Character string specifying what to plot. Default is `"qq"` which
  produces a quantile-quantile plots of the residuals.

- x:

  A vector giving the covariate values to use for residual-by- covariate
  plots (i.e., when `output = "covariate"`).

- fit:

  The fitted model from which the residuals were extracted. (Only
  required if `output = "fitted"` and `object` inherits from class
  `"resid"`.)

- distribution:

  Function that computes the quantiles for the reference distribution to
  use in the quantile-quantile plot. Default is `qnorm` which is only
  appropriate for models using a probit link function. When
  `jitter.scale = "probability"`, the reference distribution is always
  U(-0.5, 0.5). (Only required if `object` inherits from class
  `"resid"`.)

- ncol:

  Integer specifying the number of columns to use for the plot layout
  (if requesting multiple plots). Default is `NULL`.

- alpha:

  A single values in the interval \[0, 1\] controlling the opacity alpha
  of the plotted points. Only used when `nsim` \> 1.

- xlab:

  Character string giving the text to use for the x-axis label in
  residual-by-covariate plots. Default is `NULL`.

- color:

  Character string or integer specifying what color to use for the
  points in the residual vs fitted value/covariate plot. Default is
  `"black"`.

- shape:

  Integer or single character specifying a symbol to be used for
  plotting the points in the residual vs fitted value/covariate plot.

- size:

  Numeric value specifying the size to use for the points in the
  residual vs fitted value/covariate plot.

- qqpoint.color:

  Character string or integer specifying what color to use for the
  points in the quantile-quantile plot.

- qqpoint.shape:

  Integer or single character specifying a symbol to be used for
  plotting the points in the quantile-quantile plot.

- qqpoint.size:

  Numeric value specifying the size to use for the points in the
  quantile-quantile plot.

- qqline.color:

  Character string or integer specifying what color to use for the
  points in the quantile-quantile plot.

- qqline.linetype:

  Integer or character string (e.g., `"dashed"`) specifying the type of
  line to use in the quantile-quantile plot.

- qqline.size:

  Numeric value specifying the thickness of the line in the
  quantile-quantile plot.

- smooth:

  Logical indicating whether or not too add a nonparametric smooth to
  certain plots. Default is `TRUE`.

- smooth.color:

  Character string or integer specifying what color to use for the
  nonparametric smooth.

- smooth.linetype:

  Integer or character string (e.g., `"dashed"`) specifying the type of
  line to use for the nonparametric smooth.

- smooth.size:

  Numeric value specifying the thickness of the line for the
  nonparametric smooth.

- fill:

  Character string or integer specifying the color to use to fill the
  boxplots for residual-by-covariate plots when `x` is of class
  `"factor"`. Default is `NULL` which colors the boxplots according to
  the factor levels.

- resp_name:

  Character string to specify the response name that will be displayed
  in the figure.

- ...:

  Additional optional arguments to be passed onto
  [`ggplot`](https://ggplot2.tidyverse.org/reference/ggplot.html).

## Value

A `"ggplot"` object.

A `"ggplot"` object.

## Examples

``` r
# Load data
data(df1)
# Fit cumulative link model
fit <- glm(y ~ x + I(x^2), data = df1, family = binomial)
# Construct residual plots
p1 <- ggplot2::autoplot(fit, jitter.scale = "probability", output = "qq")
#> Warning: `aes_string()` was deprecated in ggplot2 3.0.0.
#> ℹ Please use tidy evaluation idioms with `aes()`.
#> ℹ See also `vignette("ggplot2-in-packages")` for more information.
#> ℹ The deprecated feature was likely used in the PAsso package.
#>   Please report the issue at <https://github.com/XiaoruiZhu/PAsso/issues>.
#> Warning: Using `size` aesthetic for lines was deprecated in ggplot2 3.4.0.
#> ℹ Please use `linewidth` instead.
#> ℹ The deprecated feature was likely used in the PAsso package.
#>   Please report the issue at <https://github.com/XiaoruiZhu/PAsso/issues>.
p2 <- ggplot2::autoplot(fit, output = "covariate", x = df1$x)
p3 <- ggplot2::autoplot(fit, output = "fitted")
```

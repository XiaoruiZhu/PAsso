# Simulated quadratic data

Data simulated from a probit model with a quadratic trend. The data are
described in Example 2 of Liu and Zhang (2017).

## Usage

``` r
data(df1)
```

## Format

A data frame with 2000 rows and 2 variables.

- `x` The predictor variable.

- `y` The response variable; an ordered factor.

## References

Liu, Dungang and Zhang, Heping. Residuals and Diagnostics for Ordinal
Regression Models: A Surrogate Approach. *Journal of the American
Statistical Association*

## Examples

``` r
head(df1)
#>   y        x
#> 1 2 5.629381
#> 2 3 6.048645
#> 3 2 3.297068
#> 4 2 3.690209
#> 5 2 5.523061
#> 6 2 4.988184
```

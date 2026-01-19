# Summary of partial association analysis

This function summarizes the partial association analysis by providing
partial association matrix, marginal association matrix, and a matrix of
the models' coefficients. The partial correlation coefficient matrix
displays the partial association between each pair of responses after
adjusting the covariates. While the marginal coefficient matrix displays
association before the adjustment.

## Usage

``` r
# S3 method for class 'PAsso'
summary(object, digits = max(3L, getOption("digits") - 2L), ...)
```

## Arguments

- object:

  A PAsso object to draw the summarized results, which includes partial
  association matrix and a matrix of the models' coefficients.

- digits:

  A default number to specify decimal digit values.

- ...:

  Additional optional arguments.

## Value

For a PAsso object, print its partial association matrix, marginal
association matrix, and a matrix of the models' coefficients.

## Examples

``` r
# See PAsso for the example.
```

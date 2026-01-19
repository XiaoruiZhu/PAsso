# generate_residuals_acat

generate_residuals_acat

## Usage

``` r
generate_residuals_acat(y, X, alpha, beta, nsim = 1)
```

## Arguments

- y:

  A vector inputs the response variable.

- X:

  A data.frame inputs the covariates.

- alpha:

  A vector provides the estimated intercepts of adjacent categories
  model. If the response has k levels, there should be k+1 numbers in
  this alpha argument with the k-1 estimated intercepts. The lower bound
  and upper bound are "-Inf" and "Inf".

- beta:

  A vector provides the estimated coefficients.

- nsim:

  A number to specify the replication of residuals.

## Value

A vector or a matrix (nsim\>1) of residuals for the adjacent categories
model.

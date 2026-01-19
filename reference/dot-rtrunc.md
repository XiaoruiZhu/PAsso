# Simulate sample from truncated distribution

A function to generate a truncated distribution. Simulate one random
sample from a standard normal distribution truncated to the left in the
middle .rtrunc(1, spec = "norm", a = -Inf, b = 0)

## Usage

``` r
.rtrunc(n, spec, a = -Inf, b = Inf, ...)
```

## Arguments

- n:

  the number of observations.

- spec:

  a character string to specify the distribution.

- a:

  lower bound.

- b:

  upper bound.

- ...:

  any other arguments that can be used for the functions of different
  distributions such as "mean", "sd" for "qnorm()".

## Value

A vector containing n random samples from the truncated distribution
"spec".

# plot3D_one

plot3D_one

## Usage

``` r
plot3D_one(
  plot_list,
  rep_SRs,
  m,
  n,
  plot_titles,
  type = c("surface3D", "contour")
)
```

## Arguments

- plot_list:

  The list to save plots. Each 3D plot is between one pair of
  covariates.

- rep_SRs:

  The surrogate responses array saved in the PAsso object.

- m:

  The index of first covariate. Since the names of covariates are empty,
  need to input index number.

- n:

  The index of second covariate.

- plot_titles:

  The title of the plot.

- type:

  A character string specifying the trace type (e.g. "surface3D",
  "contour"). "contour" creates a 2D contour plot between u and v.

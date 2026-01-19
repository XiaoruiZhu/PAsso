# 3-D P-P plot and false color level plot for the inspection of the partial association analysis

A list of 3-D P-P plots (or false color level plots when
`type = "contour"`) for the inspection of the partial association
analysis. Each plot is either 3-D P-P plot or level plot from an
empirical copula trained from the surrogate residuals of a pair of
responses.

## Usage

``` r
plot3D(object, y1, y2, ...)

# Default S3 method
plot3D(object, y1, y2, ...)

# S3 method for class 'PAsso'
plot3D(object, y1, y2, type = c("surface3D", "contour"), ...)
```

## Arguments

- object:

  A PAsso class of object.

- y1:

  A string to specify the first response for the 3D plot.

- y2:

  A string to specify the second response for the 3D plot. If either one
  of the y1 or y2 is missing. The `plot3D` will draw 3D plots for all
  pairs of responses.

- ...:

  Additional optional arguments.

- type:

  A character string specifying the trace type (e.g. "surface3D",
  "contour"). "contour" creates a 2D contour plot between u and v.

## Value

If response y1 or y2 is not specified, a list of `"plotly"` objects
includes all pairs of responses will be returned (with name "response 1
v.s. response 2" etc.). If responses y1 and y2 are specified, returns a
3D plot as `"plotly"` object.

## Details

All the plots are based on surrogate residuals generated from
`"residuals"` function in `sure`. Graphics are designed based on PAsso
and `"plotly"`.

## Examples

``` r
# Did not run this to save time
# data("ANES2016")
# PAsso_3v <- PAsso(responses = c("PreVote.num", "PID", "selfLR"),
#                   adjustments = c("income.num", "age", "edu.year"),
#                   data = ANES2016)

# plot3D(PAsso_3v, y1="PID", y2="selfLR")
# plot3D(PAsso_3v, y1="PID", y2="selfLR", type = "contour")
```

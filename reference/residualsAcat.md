# This is a function to deal with the vglm object in S4.

This is a function to deal with the vglm object in S4.

## Usage

``` r
residualsAcat(
  object,
  type = c("surrogate", "sign", "general", "deviance"),
  jitter = c("latent", "uniform"),
  jitter.uniform.scale = c("probability", "response"),
  nsim = 1L,
  ...
)
```

## Arguments

- object:

  An object of class [`vglm`](https://rdrr.io/pkg/VGAM/man/vglm.html).

- type:

  The type of residuals which should be returned. The alternatives are:
  "surrogate" (default), "sign", "general", and "deviance". Can be
  abbreviated.

  `surrogate`

  :   surrogate residuals (Liu and Zhang, 2017);

  `sign`

  :   sign-based residuals;

  `general`

  :   generalized residuals (Franses and Paap, 2001);

  `deviance`

  :   deviance residuals (-2\*loglik).

- jitter:

  A character string specifying which method to use to generate the
  surrogate response values. Current options are `"latent"` and
  `"uniform"`. Default is `"latent"`.

  `latent`

  :   latent approach;

  `uniform`

  :   jittering uniform approach.

- jitter.uniform.scale:

  A character string specifying the scale on which to perform the
  jittering whenever `jitter = "uniform"`. Current options are
  `"response"` and `"probability"`. Default is `"response"`.

- nsim:

  An integer specifying the number of replicates to use. Default is `1L`
  meaning one simulation only of residuals.

- ...:

  Additional optional arguments.

## Value

A "resid" object with attributes. It contains a vector or a matrix
(nsim\>1) of residuals for the adjacent categories model.

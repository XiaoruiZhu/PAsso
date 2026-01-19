# Extract Residuals from models of VGAM

A internal function to simulate surrogate residuals for models fitted by
[`vglm`](https://rdrr.io/pkg/VGAM/man/vglm.html) and
[`vglm`](https://rdrr.io/pkg/VGAM/man/vglm.html). Now, this one support
`"vglm"` and `"vgam"`, and adjacent categories regression model by
`"vglm"`. This one may need to update to support more models from VGAM.

## Usage

``` r
residualsVGAM(
  object,
  type = c("surrogate", "sign", "general", "deviance"),
  jitter = c("latent", "uniform"),
  jitter.uniform.scale = c("probability", "response"),
  nsim = 1L,
  ...
)
```

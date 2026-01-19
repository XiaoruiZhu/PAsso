# Changelog

## PAsso 0.1.11

### Minor improvements and fixes

1.  Change the y-axis label’s naming by removing the extra () when the
    covariate name is not specified in the diagnostic.plot() function;

2.  Allow users to specify the values alpha, shape, size in the function
    “plot”;

3.  Avoid using tidyverse as dependent package. close
    [\#11](https://github.com/XiaoruiZhu/PAsso/issues/11)

## PAsso 0.1.10

CRAN release: 2021-06-18

### Major changes

1.  Add contour plot to the function “plot3D” when argument type =
    “contour”;

### Bugs fixes

1.  Fix a warning triggered by the “.export()” in “foreach”;

## PAsso 0.1.9

CRAN release: 2021-05-07

### Minor changes

1.  Revise description of the package;

2.  Update README file;

3.  Allow users to specify the values (e.g., alpha, shape, size);

### Bugs fixes

1.  Fix Issue [\#2](https://github.com/XiaoruiZhu/PAsso/issues/2).

2.  Fix Issue [\#4](https://github.com/XiaoruiZhu/PAsso/issues/4),
    [\#5](https://github.com/XiaoruiZhu/PAsso/issues/5), and
    [\#6](https://github.com/XiaoruiZhu/PAsso/issues/6).

## PAsso 0.1.8

CRAN release: 2020-07-13

### Major changes

1.  Initialize Travis;

2.  Add cran-comments.md file to record CRAN notes;

3.  Fix errors and warnings when building PAsso locally and on Travis;

4.  Fix diagnostic.plot documentation;

5.  Clean up dependents to speed up loading of the PAsso;

### Bugs fixes

1.  Update README and other files.

## PAsso 0.1.7

### Major changes

1.  Update “PAsso()” to support adjacent categories regression model by
    using vglm(acat());

2.  Update “residuals()” to support acat in VGAM;

3.  Remove useless dataset df_AdjCat; functions “gof”;

4.  Rename attribute “boot_reps” as “draws”.

### Bugs fixes

1.  The wolfsigma method of correlation works now even when n_draws\>1.
2.  A issue when n_draws=1 in residuals.PAsso(), save array even if
    n_draws=1.

## PAsso 0.1.6

### Major changes

1.  Add ‘diagnostic.plot’ to replace ‘check.model’ for the diagnostics
    of fitted models in PAsso object. ‘autoplot’ become internal S3
    method.

2.  S3 ‘diagnostic.plot’ method can deal with classes of PAsso,
    PAsso.test, resid, clm, glm, lrm, orm, polr.

3.  Only keep one nes2016 dataset.

4.  Change ‘resids’ function to a S3 method ‘residuals’, which is
    consistent with R core ‘stats’ package.

5.  Conduct more test by ‘testthat’ to avoid bugs when methods are
    applied to different classes.

### Bug fixes

1.  Jittering on the probability scale is currently only supported for
    logit-type models.
2.  issues of ‘plot’ and ‘diagnostic.plot’.

## PAsso 0.1.5

### Major changes

1.  Replace ‘cor’ function by ‘pcaPP::cor.fk(X)’ to speed up the
    PAsso().

2.  Finalize the S3 methods ‘summary.PAsso’ and ‘print.PAsso.test’ to
    provide clean results.

### Bug fixes

1.  The digits in the print and summary are consistent.
2.  The manual is updated and S3 ‘autoplot’ is shortened.
3.  check.model for qq-plot is fixed.

## PAsso 0.1.4

### Major changes

1.  Change the name of package to “PAsso”.
2.  Change the autoplot function as a S3 method.

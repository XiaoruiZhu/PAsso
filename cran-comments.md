## Test environments
* local OS X install, R 3.6.3
* ubuntu 16.04.6 LTS (on travis-ci), R 3.6.3
* win-builder (devel and release)

## R CMD check results
There were no ERRORs, Notes, or WARNINGs. 

* [Fixed] Found the following URLs which should use \doi (with the DOI name only):
    File 'PAsso.Rd':
      https://doi.org/10.1198/jasa.2010.tm09386
      https://doi.org/10.1093/biomet/asr073
    File 'residuals.Rd':
      https://doi.org/10.1198/jasa.2010.tm09386

* [Fixed] Check: R code for possible problems, Result: NOTE
  residuals.clm: no visible global function definition for
    'residuals.vglm'

* [Fixed] checking for hidden files and directories ... NOTE
  Found the following hidden files and directories:
    .travis/.travis.yml
    .travis
  These were most likely included in error. See section 'Package
  structure' in the 'Writing R Extensions' manual.

* [Fixed] checking CRAN incoming feasibility ... NOTE
  Maintainer: 'Xiaorui (Jeremy) Zhu <zhuxiaorui1989@gmail.com>'

  
* [Fixed] checking use of S3 registration ... WARNING
  Registered S3 method from a standard package overwritten by 'PAsso':
  method        from 
  residuals.glm stats
  [An "ord" class is defined for the fitted models with ordinal response such that residuals() can recognize models!]
  
* [Fixed] Fix comment from CRAN: "You are changing the user's par() settings in your examples. Please 
reset the settings." by "dev.off()" 

* [Fixed] Add tags for for the follows:
      PAsso/man/generate_residuals_acat.Rd: \value
      PAsso/man/grid.arrange.Rd: \arguments,  \value
      PAsso/man/p_adj_cate.Rd: \value
      PAsso/man/print.Rd: \value
      PAsso/man/residualsAcat.Rd: \value
      PAsso/man/summary.Rd: \value
      PAsso/man/test.Rd: \value

## Downstream dependencies

All packages that I could install passed.

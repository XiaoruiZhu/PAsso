## Test environments
* local OS X install, R 3.6.3
* ubuntu 16.04.6 LTS (on travis-ci), R 3.6.3
* win-builder (devel and release)

## R CMD check results
There were no ERRORs or WARNINGs. 

There was 1 NOTEs:

* checking for hidden files and directories ... NOTE
  Found the following hidden files and directories:
    .travis/.travis.yml
    .travis
  These were most likely included in error. See section 'Package
  structure' in the 'Writing R Extensions' manual.

* checking CRAN incoming feasibility ... NOTE
Maintainer: 'Xiaorui(Jeremy) Zhu <zhuxiaorui1989@gmail.com>'

  New submission

  Possibly mis-spelled words in DESCRIPTION:
    Dungang (6:6)
    Irini (6:41)
    Liu (6:14, 10:27)
    Moustaki (6:47)
    Shaobo (6:19)
    Yan (6:30)
    Yu (6:34)
    Zhang (10:35)
  
* [Fixed] checking use of S3 registration ... WARNING
Registered S3 method from a standard package overwritten by 'PAsso':
  method        from 
  residuals.glm stats
  [An "ord" class is defined for the fitted models with ordinal response such that residuals() can recognize models!]
  
## Downstream dependencies

All packages that I could install passed.

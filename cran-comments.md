## Test environments
* local OS X install, R 3.6.3
* ubuntu 16.04.6 LTS (on travis-ci), R 3.6.3
* win-builder (devel and release)

## R CMD check results
There were no ERRORs or WARNINGs. 

There was 2 NOTEs:

* checking for hidden files and directories ... NOTE
  Found the following hidden files and directories:
    .travis.yml
  These were most likely included in error. See section 'Package
  structure' in the 'Writing R Extensions' manual.

* checking top-level files ... NOTE
  Non-standard files/directories found at top level:
    'cran-comments.md' 'data-raw' 


## Downstream dependencies
All packages that I could install passed.

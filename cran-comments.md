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

## Downstream dependencies
All packages that I could install passed.

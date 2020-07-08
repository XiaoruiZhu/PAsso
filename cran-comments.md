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

* Possibly mis-spelled words in DESCRIPTION:
    Dungang (6:6)
    Irini (6:41)
    Liu (6:14, 10:27)
    Moustaki (6:47)
    Shaobo (6:19)
    Yan (6:30)
    Yu (6:34)
    Zhang (10:35)
  
  Found the following (possibly) invalid URLs:
    URL: 
      From: README.md
      Message: Empty URL
    URL: https://cran.r-project.org/web/checks/check_results_PAsso.html
      From: README.md
      Status: 404
      Message: Not Found
  
  The Title field should be in title case. Current version is:
  'Assessing the Partial Association between ordinal variables'
  In title case that is:
  'Assessing the Partial Association Between Ordinal Variables'
  
  The Description field contains
    Liu and Zhang (2017) (<https://doi.org/10.1080/01621459.2017.1292915>).
  Please write DOIs as <doi:10.prefix/suffix>.

## Downstream dependencies
All packages that I could install passed.

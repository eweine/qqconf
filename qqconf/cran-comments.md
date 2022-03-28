## Test environments
* local OS X install, R 4.1.2
* win-builder (devel and release)

## R CMD check results
There were no ERRORs or WARNINGs.

There is one NOTEs, caused by re-submission to CRAN two days after a previous submission. This is necessary
to fix an issue for the package on M1 Macs.

## Downstream dependencies
I have also run R CMD check on downstream dependencies of qqconf, all packages passed


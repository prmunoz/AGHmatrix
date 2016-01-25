# AGHmatrix
`AGHmatrix` is an R package to compute A (pedigree), G (genomic-base), and H (A corrected by G) matrices for diploid and autotetraploid species to later be used in prediction models. At the moment, just the A matrix part is fully working and documented. We are working at the G and H part.

The orignal manuscript about `AGHmatrix` development and application has just been submitted to The Plant Genome journal with the following title ``AGHmatrix: R package to construct relationship matrices for autotetraploid and diploid species, a Blueberry Example``. We will keep you updated when we have news about it.

This github page is under develpment and is based upon the `OneMap` (https://github.com/augusto-garcia/onemap) software git.

# How to install

## From github

Within R, you need to install and load the package `devtools`:

```R
install.packages("devtools")
library(devtools)
```

This will allow you to automatically build and install packages from
github. If you use Windows, first install
[Rtools](https://cran.r-project.org/bin/windows/Rtools/). On a Mac,
you will need Xcode (available on the App Store). On Linux, you are
good to go.


Then, to install `AGHmatrix` from github:

```R
install_github("prmunoz/AGHmatrix")
```

Please go to the tutorial inside the vignettes folder for a complete tutorial.

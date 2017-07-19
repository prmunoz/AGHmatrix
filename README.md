# AGHmatrix

`AGHmatrix` is an [R](http://www.r-project.org) package to compute A (pedigree), G (genomic-base), and H (A corrected by G) matrices for diploid and autotetraploid species to later be used in prediction models. At the moment, just the A matrix part is fully working and documented. We are working at the G and H part.

An original manuscript about `AGHmatrix` development and application is described on [Amadeu, Rodrigo R., et al. "AGHmatrix: R Package to Construct Relationship Matrices for Autotetraploid and Diploid Species: A Blueberry Example." The Plant Genome (2016). doi:10.3835/plantgenome2016.01.0009.](https://dl.sciencesocieties.org/publications/tpg/articles/0/0/plantgenome2016.01.0009).

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

# Tutorials

You can read the _AGHmatrix_ tutorial going to the vignettes of the
installed package, or clicking below. Please, start with the overview,
that will guide you through other chapters.

[AGHmatrix Tutorial](http://htmlpreview.github.io/?https://github.com/rramadeu/aghmatrix/blob/master/inst/doc/Tutorial_AGHmatrix.html)


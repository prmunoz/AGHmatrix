# AGHmatrix

`AGHmatrix` is an [R](http://www.r-project.org) package to compute A (pedigree), G (genomic-base), and H (A corrected by G) matrices for diploid and autopolyploid species. It suports any even ploidy.

The following matrices are implemented:
### Pedigree-based relationship matrix (A matrix)

<center> 

|               | Additive                  |Non-Additive                |
|---------------|:-------------------------:|:--------------------------:|
| **Diploid**   | Henderson (1976)          |Cockerham (1954)            |
| **Polyploid** | Kerr (2012), Slater (2013)|                            |
</center>

### Molecular-based relationship matrix (G matrix) 

<center> 
  
|               | Additive                       | Non-Additive                   |
|---------------|:------------------------------:|:------------------------------:|
| **Diploid**   | Yang (2010), VanRaden (2012)   | Su (2012)                      |
| **Polyploid** | Slater (2016), VanRaden (2012) | Slater (2016), Endelman (2018) |
</center>


An original manuscript about `AGHmatrix` development and application in autotetraploids is described on [Amadeu, Rodrigo R., et al. "AGHmatrix: R Package to Construct Relationship Matrices for Autotetraploid and Diploid Species: A Blueberry Example." The Plant Genome (2016). doi:10.3835/plantgenome2016.01.0009.](https://dl.sciencesocieties.org/publications/tpg/articles/0/0/plantgenome2016.01.0009).

## How to install

### From github

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

## Tutorials

You can read the _AGHmatrix_ tutorial going to the vignettes of the
installed package, or clicking below. Please, start with the overview,
that will guide you through other chapters.

[AGHmatrix Tutorial](http://htmlpreview.github.io/?https://github.com/rramadeu/aghmatrix/blob/master/inst/doc/Tutorial_AGHmatrix.html)

## More about us
[Bluberry Breeding & Genomics Lab, University of Florida - USA](http://www.blueberrybreeding.com)

[Statistical-Genetics Lab, University of Sao Paulo - Brazil](http://www.statgen.esalq.usp.br)

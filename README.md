[![Build Status](https://travis-ci.org/rramadeu/AGHmatrix.svg?branch=master)](https://travis-ci.org/rramadeu/AGHmatrix)

# AGHmatrix

`AGHmatrix` is an [R](http://www.r-project.org) package to compute A (pedigree), G (genomic-base), and H (A corrected by G) matrices for diploid and autopolyploid species. It suports any even ploidy.

The following matrices are implemented:
### Pedigree-based relationship matrix (A matrix)

<center> 

|               | Additive                  |Non-Additive                |
|---------------|:-------------------------:|:--------------------------:|
| **Diploid**   | Henderson (1976)          |Cockerham (1954)            |
| **Polyploid** | Kerr (2012), Slater (2014)|                            |
</center>

### Molecular-based relationship matrix (G matrix) 

<center> 
  
|               | Additive                       | Non-Additive                   |
|---------------|:------------------------------:|:------------------------------:|
| **Diploid**   | Yang (2010), VanRaden (2012)   | Su (2012), Vitezica (2013)     |
| **Polyploid** | Slater (2016), VanRaden (2012) | Slater (2016), Endelman (2018) |
</center>

### Combined pedigree and molecular-based relationship matrix (H matrix)

<center> 
  
|               | Any Effect  |
|---------------|:------------------------------:|
| **Any Ploidy**   | Munoz (2014), Martini (2018)   |
</center>


An original manuscript about `AGHmatrix` development and application in autotetraploids is described on [Amadeu, R. R., C. Cellon, J. W. Olmstead, A. A. Garcia, M. F. Resende, and P. R. Muñoz, 2016 AGHmatrix: R package to construct relationship matrices for autotetraploid and diploid species: a blueberry example. The Plant Genome 9. doi:10.3835/plantgenome2016.01.0009.](https://dl.sciencesocieties.org/publications/tpg/articles/0/0/plantgenome2016.01.0009).

## How to install

### From CRAN

Within R:

```R
install.packages("AGHmatrix")
```

### From github (development version)

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

[AGHmatrix Tutorial](https://github.com/rramadeu/AGHmatrix/raw/master/inst/docs/Tutorial_AGHmatrix.pdf)

## More about us
[Bluberry Breeding & Genomics Lab, University of Florida - USA](http://www.blueberrybreeding.com)

[Statistical-Genetics Lab, University of Sao Paulo - Brazil](http://statgen.esalq.usp.br/)

## References
Amadeu, RR, et al., 2016 AGHmatrix: R package to construct relationship matrices for autotetraploid and diploid species: a blueberry example. *The Plant Genome* 9(4). https://doi.org/10.3835/plantgenome2016.01.0009

Ashraf, BH, et a., 2016 Estimating genomic heritabilities at the level of family-pool samples of perennial ryegrass using genotyping-by-sequencing. *Theoretical and Applied Genetics* 129: 45-52. https://doi.org/0.1007/s00122-015-2607-9

Endelman, JB, et al., 2018. Genetic variance partitioning and genome-wide prediction with allele dosage information in autotetraploid potato. *Genetics*, 209(1) pp. 77-87. https://doi.org/10.1534/genetics.118.300685

Hamilton, MG, et al., 2017 Computation of the inverse additive relationship matrix for autopolyploid and multiple-ploidy populations. *Theoretical and Applied Genetics*. https://doi.org/10.1007/s00122-017-3041-y

Henderson, C, 1976 A simple method for computing the inverse of a numerator relationship matrix used in prediction of breeding values. *Biometrics* pp. 69–83. https://doi.org/10.2307/2529339

Kerr, RJ, et al., 2012 Use of the numerator relation ship matrix in genetic analysis of autopolyploid species. *Theoretical and Applied Genetics* 124: 1271–1282. https://doi.org/10.1007/s00122-012-1785-y

Martini, JW, et al., 2018, The effect of the H$^{1}$ scaling factors $\tau$ and $\omega$ on the structure of H in the single-step procedure. Genetics Selection Evolution, 50(1), 16. https://doi.org/10.1186/s12711-018-0386-x

Mrode, R. A., 2014 *Linear models for the prediction of animal breeding values*. Cabi. 3rd ed.

Munoz, PR, et al., 2014 Unraveling additive from nonadditive effects using genomic relationship matrices. *Genetics* 198: 1759–1768. https://doi.org/10.1534/genetics.114.171322

R Core Team, 2016 *R*: A Language and Environment for Statistical Computing. R Foundation for Statistical Computing, Vienna, Austria.

Resende, MF, et al., 2012 Accuracy of genomic selection methods in a standard data set of loblolly pine (*Pinus taeda* l.). *Genetics* 190: 1503–1510. https://doi.org/10.1534/genetics.111.137026

Slater, AT, et al., 2014 Improving the analysis of low heritability complex traits for enhanced genetic gain in potato. *Theoretical and applied genetics* 127: 809–820. https://doi.org/10.1007/s00122-013-2258-7

Slater AT, et al., 2016 Improving genetic gain with genomic selection in autotetraploid potato. *The Plant Genome* 9.  https://doi.org/10.3835/plantgenome2016.02.0021 

Su, G, et al., 2012 Estimating additive and non-additive genetic variances and predicting genetic merits using genome-wide dense single nucleotide polymorphism markers. *PloS one* 7:e45293. https://doi.org/10.1371/journal.pone.0045293

VanRaden, P, 2008 Efficient methods to compute genomic predictions. *Journal of dairy science* 91: 4414–4423. https://doi.org/10.3168/jds.2007-0980

Vitezica, ZG, et al., 2013 On the additive and dominant variance and covariance of individuals within the genomic selection scope. *Genetics* 195: 1223–1230. https://doi.org/10.1534/genetics.113.155176

Yang, J, et al., 2010 Common snps explain a large proportion of the heritability for human height. *Nature genetics* 42: 565–569. https://doi.org/10.1038/ng.608

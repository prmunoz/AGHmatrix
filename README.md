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


An original manuscript about `AGHmatrix` development and application in autotetraploids is described on [Amadeu, R. R., C. Cellon, J. W. Olmstead, A. A. Garcia, M. F. Resende, and P. R. Muñoz, 2016 AGHmatrix: R package to construct relationship matrices for autotetraploid and diploid species: a blueberry example. The Plant Genome 9. doi:10.3835/plantgenome2016.01.0009.](https://dl.sciencesocieties.org/publications/tpg/articles/0/0/plantgenome2016.01.0009).

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

[Statistical-Genetics Lab, University of Sao Paulo - Brazil](http://statgen.esalq.usp.br/)

## References
Amadeu, R. R., C. Cellon, J. W. Olmstead, A. A. Garcia, M. F. Resende, and P. R. Muñoz, 2016 AGHmatrix: R package to construct relationship matrices for autotetraploid and diploid species: a blueberry example. The Plant Genome 9.

Cockerham (1954) as in Chapter 12 of Mrode, R. A., 2014 Linear models for the prediction of animal breeding values. Cabi. 3rd ed.

Endelman, J. B., et al., 2018. Genetic variance partitioning and genome-wide prediction with allele dosage information in autotetraploid potato. Genetics, 209(1) pp. 77-87.

Henderson, C., 1976 A simple method for computing the inverse of a numerator relationship matrix used in prediction of breeding values. Biometrics pp. 69–83.

Kerr, R. J., L. Li, B. Tier, G. W. Dutkowski, and T. A. McRae, 2012 Use of the numerator relation ship matrix in genetic analysis of autopolyploid species. Theoretical and Applied Genetics 124: 1271–1282.

Slater, A. T., G. M. Wilson, N. O. Cogan, J. W. Forster, and B. J. Hayes, 2014 Improving the analysis of low heritability complex traits for enhanced genetic gain in potato. Theoretical and applied genetics 127: 809–820.

Slater, A. T., N. O. Cogan, J. W. Forster, B. J. Hayes, and H. D. Daetwyler, 2016 Improving genetic gain with genomic selection in autotetraploid potato. The plant genome, 9(3).

Su, G., O. F. Christensen, T. Ostersen, M. Henryon, and M. S. Lund, 2012 Estimating additive and non-additive genetic variances and predicting genetic merits using genome-wide dense single nucleotide polymorphism markers. PloS one 7: e45293.

VanRaden, P., 2008 Efficient methods to compute genomic predictions. Journal of dairy science 91: 4414–4423.

Vitezica, Z. G., L. Varona, and A. Legarra, 2013 On the additive and dominant variance and covariance of individuals within the genomic selection scope. Genetics 195: 1223–1230.

Yang, J., B. Benyamin, B. P. McEvoy, S. Gordon, A. K. Henders, D. R. Nyholt, P. A. Madden, A. C. Heath, N. G. Martin, G. W. Montgomery, et al., 2010 Common snps explain a large proportion of the heritability for human height. Nature genetics 42: 565–569.

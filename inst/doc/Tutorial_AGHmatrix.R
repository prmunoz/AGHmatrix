## ----knitr_init, echo=FALSE, cache=FALSE---------------------------------
library(knitr)
library(rmarkdown)

knitr::opts_chunk$set(collapse = TRUE,
                      comment = "#>",
                      fig.width = 6,
                      fig.height = 6,
                      fig.align = "center",
                      dev = "png",
                      dpi = 36,
                      cache = TRUE)


## ---- echo=FALSE, results='hide'-----------------------------------------
library(AGHmatrix)

## ---- eval=FALSE---------------------------------------------------------
#  ## Installing and loading devtools
#  install.packages("devtools")
#  library(devtools)
#  
#  ## Installing from github the fullsibQTL
#  install_github("prmunoz/AGHmatrix")

## ---- eval=TRUE----------------------------------------------------------
library(AGHmatrix)

## ------------------------------------------------------------------------
data(ped.mrode)
ped.mrode

## ------------------------------------------------------------------------
Amatrix(ped.mrode, ploidy=2)

## ---- eval=FALSE---------------------------------------------------------
#  #Computing non-additive relationship matrix considering diploidy:
#  Amatrix(ped.mrode, ploidy=2, dominance=TRUE)
#  
#  #Computing additive relationship matrix considering autotetraploidy:
#  Amatrix(ped.mrode, ploidy=4)
#  
#  #Computing additive relationship matrix considering autooctaploidy:
#  Amatrix(ped.mrode, ploidy=8)
#  
#  #Computing additive relationship matrix considering autotetraploidy and double-reduction of 10%:
#  Amatrix(ped.mrode, ploidy=4, w=0.1)
#  
#  #Computing additive relationship matrix considering autohexaploidy and double-reduction of 10%:
#  Amatrix(ped.mrode, ploidy=6, w=0.1)

## ---- eval=FALSE---------------------------------------------------------
#  ?Amatrix

## ------------------------------------------------------------------------
data(snp.pine)
str(snp.pine)

## ---- eval=FALSE---------------------------------------------------------
#  #Computing the additive relationship matrix based on VanRaden 2008
#  G_VanRaden <- Gmatrix(SNPmatrix=snp.pine, missingValue=-9, maf=0.05, method="VanRaden")
#  
#  #Computing the additive relationship matrix based on Yang 2010
#  G_Yang <- Gmatrix(SNPmatrix=snp.pine, missingValue=-9, maf=0.05, method="Yang")
#  
#  #Computing the dominance relationship matrix based on Su 2012
#  G_Su <- Gmatrix(SNPmatrix=snp.pine, missingValue=-9, maf=0.05, method="Su")
#  
#  #Computing the dominance relationship matrix based on Vitezica 2013
#  G_Vitezica <- Gmatrix(SNPmatrix=snp.pine, missingValue=-9, maf=0.05, method="Vitezica")

## ---- eval=FALSE---------------------------------------------------------
#  ?Gmatrix

## ------------------------------------------------------------------------
inds <- 10
markers <- 100
markersdata <- matrix(sample(x=0:4, size=inds*markers, replace=TRUE), nrow=inds, ncol=markers)

## ---- eval=FALSE---------------------------------------------------------
#  #Computing the additive relationship matrix based on VanRaden 2008 and adapted by Ashraf 2016
#  G_VanRaden <- Gmatrix(markersdata, method="VanRaden", ploidy=4)
#  
#  #Computing the dominance (digenic) matrix based on Endelman 2018 (Eq. 19)
#  G_Digenic <- Gmatrix(markersdata, method="Endelman", ploidy=4)
#  
#  #Computing the full-autopolyploid matrix based on Slater 2016 (Eq. 8 and 9)
#  G_FullAutopolyploid <- Gmatrix(markersdata, method="Slater", ploidy=4)
#  
#  #Computing the pseudodiploid matrix based on Slater 2016 (Eq. 5, 6, and 7)
#  G_Pseudodiploid <- Gmatrix(markersdata, method="VanRaden", ploidy=4, pseudo.diploid=TRUE)

## ---- eval=FALSE---------------------------------------------------------
#  ?Gmatrix

## ---- eval=FALSE---------------------------------------------------------
#  #Computing the additive relationship matrix based on Yang 2010
#  Gmat <- Gmatrix(SNPmatrix=snp.pine, missingValue=-9, maf=0.05, method="VanRaden")
#  
#  #Setting an A matrix structure
#  Amat <- diag(926) #we have 926 individuals
#  
#  #A matrix has to have the same row and col names of G
#  rownames(Amat) <- colnames(Amat) <- rownames(Gmat)
#  
#  #Computing H matrix
#  Hmat <- Hmatrix(A=Amat, G=Gmat, markers=snp.pine, missingValue=-9, maf=0.05)

## ---- eval=FALSE---------------------------------------------------------
#  A <- Gmatrix(SNPmatrix=snp.pine, method="VanRaden", missingValue=-9, maf=0.05)
#  D <- Gmatrix(SNPmatrix=snp.pine, method="Vitezica", missingValue=-9,maf=0.05)

## ---- eval=FALSE---------------------------------------------------------
#  #Additive-by-Additive Interactions
#  A_A <- A*A
#  #Dominance-by-Additive Interactions
#  D_A <- D*A
#  #Dominance-by-Dominance Interactions
#  D_D <- D*D

## ---- eval=FALSE---------------------------------------------------------
#  #Additive-by-Additive-by-Additive Interactions
#  A_A_A <- A*A*A
#  #Additive-by-Additive-by-Dominance Interactions
#  A_A_D <- A*A*D
#  #Additive-by-Dominance-by-Dominance Interactions
#  A_D_D <- A*D*D
#  #Dominance-by-Dominance-by-Dominance Interactions
#  D_D_D <- D*D*D

## ---- eval=FALSE---------------------------------------------------------
#  #Computing the matrix
#  A <- Amatrix(data=ped.mrode, ploidy=4, w=0.1)
#  
#  #Building its unique inverse
#  Ainv <- solve(A)

## ---- eval=FALSE---------------------------------------------------------
#  #Computing the matrix
#  Gmat <- Gmatrix(SNPmatrix=snp.pine, missingValue=-9, maf=0.05, method="VanRaden")
#  
#  #Building its Moore-Penrose generalized inverse
#  Ginv <- MASS::ginv(Gmat)

## ---- eval=FALSE---------------------------------------------------------
#  #Loading the data example
#  data(ped.mrode)
#  
#  #Computing the matrix
#  A <- Amatrix(data=ped.mrode, ploidy=4, w=0.1)
#  
#  #Building its inverse
#  Ainv <- solve(A)
#  
#  #Exporting it. The function "formatmatrix" will convert it and save in your working directory
#  formatmatrix(Ainv, round.by=12, exclude.0=TRUE, name="Ainv")

## ---- eval=FALSE---------------------------------------------------------
#  #Loading the data example
#  data(ped.mrode)
#  
#  #Determining a double reduction range
#  double.red<-seq(0,0.2,0.05)
#  
#  #Extracting the length of double.red
#  n<-length(double.red)
#  
#  #Making the loop
#  for(i in 1:n){
#    A<-Amatrix(data=ped.mrode,
#               ploidy=4,
#               w=double.red[i])
#    #Computing the inverse
#    A.inv<-solve(A)
#  
#    #Exporting as csv
#    formatmatrix(data=A.inv,
#                 name=paste("Ainv_",double.red[i],sep=""),
#                 round.by=12,
#                 exclude.0=TRUE)
#  }
#  

## ------------------------------------------------------------------------
sessionInfo()


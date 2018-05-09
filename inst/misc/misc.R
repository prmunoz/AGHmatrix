library(AGHmatrix)
Gmatrix

inds <- 10
markers <- 100
markersdata <- matrix(sample(x=0:4, size=inds*markers, replace=TRUE), nrow=inds, ncol=markers)

markersdata
G1 <-Gmatrix(markersdata,ploidy=4,method="VanRaden")
G2 <-Gmatrix(markersdata,ploidy=4,method="Endelman")

SNPmatrix<-markersdata

inds <- 3
markers <- 4
markersdata <- matrix(sample(x=0:4, size=inds*markers, replace=TRUE), nrow=inds, ncol=markers)

Frequency <- colSums(SNPmatrix)/(nrow(SNPmatrix)*ploidy)
Frequency <- cbind(Frequency,1-Frequency)

Frequency <- cbind(apply(SNPmatrix, 2, function(x) alelleFreq(x,0))
                   , apply(SNPmatrix, 2, function(x) alelleFreq(x, 2)))

Q1 <- matrix(rep(c(9,4,25,16),3),nrow=3,byrow=TRUE)/144*6
Q2 <- matrix(c(0,4,10,4,
               9,0,0,12,
               0,0,15,0),nrow=3,byrow=TRUE)*3/12
Q3 <- matrix(c(0,2,2,0,
               6,0,0,6,
               0,0,6,0),nrow=3,byrow=TRUE)*0.5
Q<-Q1-Q2+Q3
tcrossprod(Q)/drop((6*crossprod(c(3/12,2/12,5/12,4/12)^2,c(1-c(3/12,2/12,5/12,4/12))^2)))

M1 <- matrix(c(0,2,2,1,
              3,0,0,3,
              0,0,3,0),nrow=3,byrow=TRUE)

M2 <-matrix(c(4,2,2,3,
             1,4,4,1,
             4,4,1,4),nrow=3,byrow=TRUE)

Gmatrix(M1,ploidy=4)
Gmatrix(M2,ploidy=4)

## Including Data from Endelman 2018




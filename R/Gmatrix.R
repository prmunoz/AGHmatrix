#########################################
# 									
# Package: AGHmatrix 							
# 									
# File: Gmatrix.R
# Contains: Gmatrix slater_par						
# 									
# Written by Rodrigo Rampazo Amadeu 			
# Collaborators: Marcio Resende Jr and Leticia AC Lara
# 									
# First version: Feb-2014 					
# Last update: 22-Feb-2017 						
# License: GPL-3	
# 									
#########################################

#' Construction of Relationship Matrix G
#'
#' Given a matrix (individual x markers), a method, a missing value, and a maf threshold, return a additive or dominance relationship matrix. Uses pseudodiploid model (Slater, 2016) on polyploid data for methods "Yang", "VanRaden", "Su" and "Vitezica".
#'
#' @param SNPmatrix matrix (n x m), where n is is individual names and m is marker names (coded inside the matrix as 0, 1, 2, missingValue). 
#' @param method "Yang" or "VanRaden" for marker-based additive relationship matrix. "Su" or "Vitezica" for marker-based dominance relationship matrix. "Slater" for full-autopolyploid model including non-additive effects. "MarkersMatrix" for a matrix with the amount of shared markers between individuals (3). Default="Yang".
#' @param missingValue missing value in data. Default=-9.
#' @param maf max of missing data accepted to each marker. Default=0.05.
#' @param verify.posdef verify if the resulting matrix is positive-definite. Default=TRUE.
#' @param ploidy data ploidy (an even number between 2 and 20). Default=2.
#'
#' @return Matrix with the marker-bases relationships between the individuals
#'
#' @examples
#' data(snp.pine)
#' #Verifying if data is coded as 0,1,2 and missing value.
#' str(snp.pine)
#' #Build Gmatrix
#' Gmatrix.Yang <- Gmatrix(snp.pine,method="Yang",missingValue=-9,maf=0.05) 
#' 
#' @author Rodrigo R Amadeu \email{rramadeu@@gmail.com}, Marcio Resende Jr, and Letícia Aparecida de Castro Lara
#' 
#' @references \emph{Slater, A.T., Cogan, N.O., Forster, J.W., Hayes, B.J., Daetwyler, H.D, 2016. Improving genetic gain with genomic selection in autotetraploid potato. The Plant Genome 9(3), pp.1-15.}
#' @references \emph{Su, G., Christensen, O.F., Ostersen, T., Henryon, M. and Lund, M.S., 2012. Estimating additive and non-additive genetic variances and predicting genetic merits using genome-wide dense single nucleotide polymorphism markers. PloS one, 7(9), p.e45293.}
#' @references \emph{VanRaden, P.M., 2008. Efficient methods to compute genomic predictions. Journal of dairy science, 91(11), pp.4414-4423.}
#' @references \emph{Vitezica, Z.G., Varona, L. and Legarra, A., 2013. On the additive and dominant variance and covariance of individuals within the genomic selection scope. Genetics, 195(4), pp.1223-1230.}
#' @references \emph{Yang, J., Benyamin, B., McEvoy, B.P., Gordon, S., Henders, A.K., Nyholt, D.R., Madden, P.A., Heath, A.C., Martin, N.G., Montgomery, G.W. and Goddard, M.E., 2010. Common SNPs explain a large proportion of the heritability for human height. Nature genetics, 42(7), pp.565-569.}
#' 
#' @export

Gmatrix <- function (SNPmatrix = NULL, method = "Yang", missingValue = -9, maf = 0, verify.posdef = TRUE, ploidy=2){
  Time = proc.time()
  
  if (!is.na(missingValue)) {
    m <- match(SNPmatrix, missingValue, 0)
    SNPmatrix[m > 0] <- NA
  }
  
  check_Gmatrix_data(SNPmatrix=SNPmatrix,method=method,ploidy=ploidy)
  
  NumberMarkers <- ncol(SNPmatrix)
  nindTotal <- colSums(!is.na(SNPmatrix))
  nindAbs <- max(nindTotal)
  cat("\tNumber of Individuals:", max(nindTotal), "\n")
  cat("\tNumber of Markers:", NumberMarkers, "\n")
  
  if(method=="Slater"){
    P <- colSums(SNPmatrix,na.rm = TRUE)/nrow(SNPmatrix)
    SNPmatrix[,which(P>ploidy/2)] <- ploidy-SNPmatrix[,which(P>(ploidy/2))]
    SNPmatrix <- slater_par(SNPmatrix,ploidy=ploidy)
    NumberMarkers <- ncol(SNPmatrix)
    Frequency <- colSums(SNPmatrix,na.rm=TRUE)/nrow(SNPmatrix)
    FreqP <- matrix(rep(Frequency, each = nrow(SNPmatrix)), 
                    ncol = ncol(SNPmatrix))
  }
  
  if(ploidy==2){
    alelleFreq <- function(x, y) {
      (2 * length(which(x == y)) + length(which(x == 1)))/(2 * 
                                                             length(which(!is.na(x))))
    }
    Frequency <- cbind(apply(SNPmatrix, 2, function(x) alelleFreq(x,0))
                       , apply(SNPmatrix, 2, function(x) alelleFreq(x, 2)))
    
    if (any(Frequency[, 1] <= maf) & maf != 0) {
      cat("\t", length(which(Frequency[, 1] <= maf)), "markers dropped due to maf cutoff of", maf, "\n")
      SNPmatrix <- SNPmatrix[,-which(Frequency[, 1] <= maf)]
      cat("\t", ncol(SNPmatrix), "markers kept \n")
      Frequency <- as.matrix(Frequency[-which(Frequency[,1] <= 
                                                maf), ])
      NumberMarkers <- ncol(SNPmatrix)
    }
    FreqP <- matrix(rep(Frequency[, 2], each = nrow(SNPmatrix)), 
                    ncol = ncol(SNPmatrix))
  }
  
  if(ploidy>2 && method!="Slater"){## Uses Pseudodiploid model
    P <- colSums(SNPmatrix,na.rm = TRUE)/nrow(SNPmatrix)
    SNPmatrix[,which(P>ploidy/2)] <- ploidy-SNPmatrix[,which(P>(ploidy/2))]
    Frequency <- colSums(SNPmatrix,na.rm=TRUE)/(ploidy*nrow(SNPmatrix))
    Frequency <- cbind(1-Frequency,Frequency)
    FreqP <- matrix(rep(Frequency[, 2], each = nrow(SNPmatrix)), 
                    ncol = ncol(SNPmatrix))
    SNPmatrix[SNPmatrix %in% c(1:(ploidy-1))] <- 1
    SNPmatrix[SNPmatrix==ploidy] <- 2
  }
  
  if (method == "MarkersMatrix") {
    Gmatrix <- !is.na(SNPmatrix)
    Gmatrix <- tcrossprod(Gmatrix, Gmatrix)
    return(Gmatrix)
  }
  if (method == "VanRaden") {
    TwoPQ <- 2 * t(Frequency[, 1]) %*% Frequency[, 2]
    SNPmatrix <- SNPmatrix- 2 * FreqP
    SNPmatrix[is.na(SNPmatrix)] <- 0
    Gmatrix <- (tcrossprod(SNPmatrix, SNPmatrix))/as.numeric(TwoPQ)
  }
  if (method == "Yang") {
    FreqPQ <- matrix(rep(2 * Frequency[, 1] * Frequency[, 
                                                        2], each = nrow(SNPmatrix)), ncol = ncol(SNPmatrix))
    G.all <- (SNPmatrix^2 - (1 + 2 * FreqP) * SNPmatrix + 
                2 * (FreqP^2))/FreqPQ
    G.ii <- as.matrix(colSums(t(G.all), na.rm = T))
    SNPmatrix <- (SNPmatrix - (2 * FreqP))/sqrt(FreqPQ)
    G.ii.hat <- 1 + (G.ii)/NumberMarkers
    SNPmatrix[is.na(SNPmatrix)] <- 0
    Gmatrix <- (tcrossprod(SNPmatrix, SNPmatrix))/NumberMarkers
    diag(Gmatrix) <- G.ii.hat
  }
  if (method == "Su"){
    TwoPQ <- 2*(FreqP)*(1-FreqP)
    SNPmatrix[SNPmatrix==2 | SNPmatrix==0] <- 0
    SNPmatrix <- SNPmatrix - TwoPQ
    SNPmatrix[is.na(SNPmatrix)] <- 0
    Gmatrix <- tcrossprod(SNPmatrix,SNPmatrix)/
      sum(TwoPQ[1,]*(1-TwoPQ[1,]))        
  }
  if (method == "Vitezica"){
    TwoPQ <- 2*(FreqP[1,])*(1-FreqP[1,])
    SNPmatrix[is.na(SNPmatrix)] <- 0
    SNPmatrix <- (SNPmatrix==0)*-2*(FreqP^2) +
      (SNPmatrix==1)*2*(FreqP)*(1-FreqP) +
      (SNPmatrix==2)*-2*((1-FreqP)^2)
    Gmatrix <- tcrossprod(SNPmatrix,SNPmatrix)/sum(TwoPQ^2)
  }
  
  if (method == "Slater"){
    drop.alleles <- which(Frequency==0)
    Frequency <- Frequency[-drop.alleles]
    SNPmatrix <- SNPmatrix[,-drop.alleles]
    FreqP <- FreqP[,-drop.alleles]
    FreqPQ <- matrix(rep(Frequency * (1-Frequency),
                         each = nrow(SNPmatrix)),
                     ncol = ncol(SNPmatrix))
    SNPmatrix[which(is.na(SNPmatrix))] <- 0
    G.ii <- (SNPmatrix^2 - (2 * FreqP) * SNPmatrix + FreqP^2)/FreqPQ
    G.ii <- as.matrix(colSums(t(G.ii), na.rm = T))
    G.ii <- 1 + (G.ii)/NumberMarkers
    SNPmatrix <- (SNPmatrix - (FreqP))/sqrt(FreqPQ)
    SNPmatrix[is.na(SNPmatrix)] <- 0
    Gmatrix <- (tcrossprod(SNPmatrix, SNPmatrix))/NumberMarkers
    diag(Gmatrix) <- G.ii
  }
  if (verify.posdef) {
    e.values <- eigen(Gmatrix, symmetric = TRUE)$values
    indicator <- length(which(e.values <= 0))
    if (indicator > 0) 
      cat("\t Matrix is NOT positive definite. It has ", indicator, 
          " eigenvalues <= 0 \n \n")
  }
  Time = as.matrix(proc.time() - Time)
  cat("Completed! Time =", Time[3], " seconds \n")
  gc()
  return(Gmatrix)
}

## Internal Functions ##

# Coding SNPmatrix as Slater (2016) Full autotetraploid model including non-additive effects (Presence/Absence per Genotype per Marker)
slater_par <- function(X,ploidy){
  prime.index <- c(3,5,7,11,13,17,19,23,29,31,37,
                   41,43,47,53,59,61,67,71,73,79)
  
  NumberMarkers <- ncol(X)
  nindTotal <- nrow(X)
  X <- X+1
  
  ## Breaking intervals to use less RAM
  temp <- seq(1,NumberMarkers,10000)
  temp <- cbind(temp,temp+9999)
  temp[length(temp)] <- NumberMarkers
  prime.index <- prime.index[1:(ploidy+1)]
  
  ## Uses Diagonal (which is Sparse mode, uses less memmory)
  for(i in 1:nrow(temp)){
    X.temp <- X[,c(temp[i,1]:temp[i,2])]
    NumberMarkers <- ncol(X.temp)
    X.temp <- X.temp %*% t(kronecker(diag(NumberMarkers),prime.index))
    X.temp[which(as.vector(X.temp) %in%
                   c(prime.index*c(1:(ploidy+1))))] <- 1
    X.temp[X.temp!=1] <- 0
    if(i==1){
      X_out <- X.temp
    }else{
      X_out <- cbind(X_out,X.temp)
    }   
  }
  gc()
  return(X_out)
}

# Internal function to check input Gmatrix arguments
check_Gmatrix_data <- function(SNPmatrix,ploidy,method){
  if (is.null(SNPmatrix)) {
    stop(deparse("Please define the variable SNPdata"))
  }
  if (all(method != c("Yang", "VanRaden", "Slater", "Su", "Vitezica", "MarkersMatrix"))) {
    stop("Method to build Gmatrix has to be either `Yang` or `VanRaden` for marker-based additive relationship matrix, or `Su` or `Vitezica` for marker-based dominance relationship matrx, or `MarkersMatrix` for matrix with amount of shared-marks by individuals pairs")
  }
  if(class(SNPmatrix)!="matrix"){
    cat("SNPmatrix class is:",class(SNPmatrix),"\n")
    stop("SNPmatrix class must be matrix. Please verify it.")
  }
  
  if( ploidy > 20 | (ploidy %% 2) != 0)
    stop(deparse("Only even ploidy from 2 to 20"))
  
  t <- max(SNPmatrix,na.rm = TRUE)
  if( t > ploidy )
    stop(deparse("Check your data, it has values above ploidy number"))
  
  t <- min(SNPmatrix,na.rm=TRUE)
  if( t < 0 )
    stop(deparse("Check your data, it has values under 0"))
  
  if(prod(SNPmatrix == round(SNPmatrix),na.rm = TRUE)==0)
    stop(deparse("Check your data, it has not integer values"))
}
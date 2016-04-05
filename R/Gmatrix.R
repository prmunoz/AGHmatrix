#########################################################################
# 									#
# Package: AGHmatrix 							#
# 									#
# File: Gmatrix.R 							#
# Contains: Gmatrix 							#
# 									#
# Written by Rodrigo Rampazo Amadeu 					#
# 									#
# First version: Feb-2014 						#
# Last update: 31-Mar-2016 						#
# License: GNU General Public License version 2 (June, 1991) or later 	#
# 									#
#########################################################################

#' Construction of Relationship Matrix G
#'
#' Given a matrix (markers x individual) choosing a method ((1) Powell or (2) VanRaden), a missing value, and a maf threshold, return a G matrix. Also the function can calculate the number of markers shared two-by-two individuals (method 3).
#'
#' @param SNPdata matrix (n x m), where n is marker information (coded as 0,1,2,NA) and m is individual names (a string). 
#' @param method "Powell" (1), "VanRaden" (2), shared markers between individual (3). Default=2.
#' @param missingValue missing value in data. Default=NA.
#' @param maf max of missing data accepted to each marker. Default=0.
#' @param verify.posdef verify if the resulting matrix is positive-definite. Default=TRUE.
#'
#' @return Matrix with the Relationship between the individuals
#'
#' @examples
#' data(snp.table)
#' #Verifying if data is coded as 0,1,2 and missing value.
#' snp.table
#' #Build Gmatrix
#' Gmatrix(snp.table,missingValue=-9)
#'
#' @author Rodrigo R Amadeu, \email{rramadeu@@gmail.com}
#'
#' @references \emph{Van Raden, P. M. (2008). Efficient methods to compute genomic predictions. Journal of Dairy Science, 91(11), 4414-4423.}
#' @references \emph{Powell, J. E. (2010). Reconciling the analysis of IBD and IBS in complex trait studies. Nature Reviews Genetics, 11(11), 800-805.}
#'
#' @export

Gmatrix = function( SNPdata = NULL,
                    method = 2,
                    missingValue = NA,
                    maf = 0,
                    verify.posdef=TRUE){

    Time = proc.time()

    if(is.null(SNPdata)){
        stop(deparse("Please define the variable SNPdata"))
    }
    if (missing(missingValue)){ stop(deparse("missingValue not defined"))}
    
    if(method==1)
        method="Powell"
    if(method==2)
        method="VanRaden"
    if(method==3)
        method="MarkersMatrix"

    if (all(method != c("Powell","VanRaden","MarkersMatrix"))){ stop("Method to build Gmatrix has to be either (1) Powell or (2) VanRaden or (3) for the Amount of Markers Between Individuals  Matrix") }

    SNPmatrix = as.matrix(SNPdata)

    if (!is.na(missingValue)){
        m <- match(SNPmatrix, missingValue ,0)
        SNPmatrix[m > 0] <- NA
    }

    NumberMarkers = nrow(SNPmatrix)
    nindTotal = rowSums(!is.na(SNPmatrix))
    nindAbs = max(nindTotal)

    cat("	Number of Individuals:", max(nindTotal),"\n\n")
    cat("	Number of Markers:", NumberMarkers ,"\n")

    alelleFreq = rowSums(SNPmatrix,na.rm=TRUE)/(2*nindTotal)

    ##############################################
    #Swap major with minor alelle (0,1,2 to 2,1,0)
    major.minor <- which(alelleFreq < 0.5)
    if( length(major.minor>0) )
        SNPmatrix[major.minor,] <- 1+ (SNPmatrix[major.minor,] - 1)*-1
    ##############################################
   
    if( maf == 0){
        exclude <- c(which(alelleFreq==0))

        if(length(exclude)>0){
            SNPmatrix <- SNPmatrix[-exclude,]
            alelleFreq <- alelleFreq[-exclude]
        }
    }
    
    if( maf > 0){
        exclude <- c(which(alelleFreq<maf),which(alelleFreq>(1-maf)))
        if(length(exclude)>0){
            SNPmatrix <- SNPmatrix[-exclude,]
            alelleFreq = alelleFreq[-exclude]
            NumberMarkers = nrow(SNPmatrix)
        }
        cat("	Number of Markers after maf:", NumberMarkers ,"\n")
    }

    if(method=="MarkersMatrix"){ #calculate the number of markers for each ind pair
        Gmatrix <- !is.na(SNPmatrix)
        Gmatrix <- crossprod(Gmatrix,Gmatrix)
        return(Gmatrix)
    }

    NumberMarkers = nrow(SNPmatrix)

    Frequency = cbind((1-alelleFreq),alelleFreq)
    FreqP<-alelleFreq

#    Pmatrix <- matrix(rep(alelleFreq,max(nindTotal)),nrow=max(nindTotal),byrow=TRUE)
    Pmatrix <- matrix(rep(alelleFreq,max(nindTotal)),ncol=max(nindTotal))
    

    if( method=="VanRaden"){
        Zmatrix <- SNPmatrix-2*Pmatrix
        Zmatrix[is.na(Zmatrix)] = 0 #Missing Values
        Gmatrix <- crossprod(Zmatrix,Zmatrix)
        sumpi <- sum(alelleFreq*(1-alelleFreq))
        Gmatrix <- Gmatrix/(2*sumpi)
    }

    if( method == "Powell"){
        SNP.Visscher <- SNPmatrix
        FreqPQ = matrix(rep(2*Frequency[,1] * Frequency[,2],nindAbs),ncol=nindAbs)
        G.all = (SNP.Visscher^2 - (1+2*Pmatrix)*SNP.Visscher + 2*Pmatrix^2)/FreqPQ
        SNP.Visscher= (SNP.Visscher-(2*Pmatrix))/sqrt(FreqPQ)
        SNP.Visscher[is.na(SNP.Visscher)] <- 0
        Gmatrix = (crossprod(SNP.Visscher,SNP.Visscher))/max(NumberMarkers)
      
        G.ii = as.matrix(colSums(G.all,na.rm = T))
        G.ii.hat = 1+(G.ii)/max(NumberMarkers)
        diag(Gmatrix) = G.ii.hat
    }
    
   ##############################################
   #Verifying if it is positive definite
    if( verify.posdef ){
        e.values<-eigen(Gmatrix,symmetric=TRUE)$values
        indicator<-length(which(e.values <= 0))
        if( indicator > 0)
            cat("Matrix is NOT posiive definite. It has ", indicator," eigenvalues <= 0")
    }
     ##############################################

    Time = as.matrix(proc.time()-Time)
    cat("Completed! Time =", Time[3]," seconds \n")
    return(Gmatrix)
}

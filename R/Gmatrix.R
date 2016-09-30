#########################################################################
# 									#
# Package: AGHmatrix 							#
# 									#
# File: Gmatrix.R 							#
# Contains: Gmatrix 							#
# 									#
# Written by Rodrigo Rampazo Amadeu and Marcio Resende Jr.	#
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
#' @author Rodrigo R Amadeu, \email{rramadeu@@gmail.com} and Marcio Resende Jr.
#'
#' @references \emph{Van Raden, P. M. (2008). Efficient methods to compute genomic predictions. Journal of Dairy Science, 91(11), 4414-4423.}
#' @references \emph{Powell, J. E. (2010). Reconciling the analysis of IBD and IBS in complex trait studies. Nature Reviews Genetics, 11(11), 800-805.}
#'
#' @export

Gmatrix <- function (SNPdata = NULL, method = 2, missingValue = NA, maf = 0, verify.posdef = TRUE) 
{
    Time = proc.time()
    if (is.null(SNPdata)) {
        stop(deparse("Please define the variable SNPdata"))
    }
    if (missing(missingValue)) {
        stop(deparse("missingValue not defined"))
    }
    if (method == 1) 
        method = "Powell"
    if (method == 2) 
        method = "VanRaden"
    if (method == 3) 
        method = "MarkersMatrix"
    if (all(method != c("Powell", "VanRaden", "MarkersMatrix"))) {
        stop("Method to build Gmatrix has to be either (1) Powell or (2) VanRaden or (3) for the Amount of Markers Between Individuals  Matrix")
    }
 SNPmatrix = t(matrix(as.numeric(SNPdata),ncol = ncol(SNPdata),dimnames = list(c(rownames(SNPdata)),c(colnames(SNPdata)))))
    if (!is.na(missingValue)) {
        m <- match(SNPmatrix, missingValue, 0)
        SNPmatrix[m > 0] <- NA
    }
    NumberMarkers = nrow(SNPmatrix)
    nindTotal = rowSums(!is.na(SNPmatrix))
    nindAbs = max(nindTotal)
    cat("\tNumber of Individuals:", max(nindTotal), "\n\n")
    cat("\tNumber of Markers:", NumberMarkers, "\n")
    alelleFreq = function(x,y){(2*length(which(x==y))+length(which(x==1)))/(2*length(which(!is.na(x))))}
    
    Frequency = cbind(apply(SNPmatrix,1,function(x) alelleFreq(x,0)),apply(SNPmatrix,1,function(x) alelleFreq(x,2)))
    if (any(Frequency[,1]<= maf) & maf != 0){
        SNPmatrix = SNPmatrix[-which(Frequency[,1]<= maf),]
        Frequency = as.matrix(Frequency[-which(Frequency[,1]<= maf),])
        NumberMarkers = nrow(SNPmatrix)
    }	
    Sum2pq = 2*t(Frequency[,1])%*%Frequency[,2]
    FreqP = matrix(rep(Frequency[,2],each = ncol(SNPmatrix)),ncol = nrow(SNPmatrix))
    if (method == "MarkersMatrix") {
        Gmatrix <- !is.na(SNPmatrix)
        Gmatrix <- crossprod(Gmatrix, Gmatrix)
        return(Gmatrix)
    }
    if (method == "VanRaden") {
	SNP_VanRaden = as.matrix(t(SNPmatrix))
        SNP_VanRaden = SNP_VanRaden - 2*FreqP
        SNP_VanRaden[is.na(SNP_VanRaden)] <- 0
        Gmatrix = (crossprod(t(SNP_VanRaden),t(SNP_VanRaden)))/as.numeric(Sum2pq)
    }
    if (method == "Powell") {
        
		SNP_Visscher = as.matrix(t(SNPmatrix))		
		FreqPQ = matrix(rep(2*Frequency[,1] * Frequency[,2],each = ncol(SNPmatrix)),ncol = nrow(SNPmatrix))			
		G_all = (SNP_Visscher^2 - (1+2*FreqP)*SNP_Visscher + 2*FreqP^2)/FreqPQ
		G_ii = as.matrix(colSums(t(G_all),na.rm = T))
		SNP_Visscher= (SNP_Visscher-(2*FreqP))/sqrt(FreqPQ)
		G_ii_hat = 1+(G_ii)/NumberMarkers
		SNP_Visscher[is.na(SNP_Visscher)] <- 0
		Gmatrix = (crossprod(t(SNP_Visscher),t(SNP_Visscher)))/NumberMarkers
		diag(Gmatrix) = G_ii_hat
    }
    if (verify.posdef) {
        e.values <- eigen(Gmatrix, symmetric = TRUE)$values
        indicator <- length(which(e.values <= 0))
        if (indicator > 0) 
            cat("Matrix is NOT posiive definite. It has ", indicator, 
                " eigenvalues <= 0")
    }
    Time = as.matrix(proc.time() - Time)
    cat("Completed! Time =", Time[3], " seconds \n")
    return(Gmatrix)
}

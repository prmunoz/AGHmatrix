#########################################################################
# 									#
# Package: AGHmatrix 							#
# 									#
# File: Gmatrix.R 							#
# Contains: Gmatrix 							#
# 									#
# Written by Rodrigo Rampazo Amadeu and Marcio Resende Jr.		#
# 									#
# First version: Feb-2014 						#
# Last update: 06-Oct-2016 						#
# License: GNU General Public License version 2 (June, 1991) or later 	#
# 									#
#########################################################################

#' Construction of Relationship Matrix G
#'
#' Given a matrix (individual x markers), a method ((1) Powell or (2) VanRaden), a missing value, and a maf threshold, return a G matrix. Also the function can calculate the number of markers shared two-by-two individuals (method 3).
#'
#' @param SNPmatrix matrix (n x m), where n is is individual names and m is marker names (coded inside the matrix as 0,1,2,missingValue). 
#' @param method "Powell" (1), "VanRaden" (2), shared markers between individual (3). Default=2.
#' @param missingValue missing value in data. Default=-9.
#' @param maf max of missing data accepted to each marker. Default=0.
#' @param verify.posdef verify if the resulting matrix is positive-definite. Default=TRUE.
#'
#' @return Matrix with the Relationship between the individuals
#'
#' @examples
#' data(snp.pine)
#' #Verifying if data is coded as 0,1,2 and missing value.
#' str(snp.pine)
#' #Build Gmatrix
#' PineGmatrix <- Gmatrix(snp.pine,missingValue=-9,maf=0.05)
#'
#' @author Rodrigo R Amadeu, \email{rramadeu@@gmail.com} and Marcio Resende Jr.
#'
#' @references \emph{Van Raden, P. M. (2008). Efficient methods to compute genomic predictions. Journal of Dairy Science, 91(11), 4414-4423.}
#' @references \emph{Powell, J. E. (2010). Reconciling the analysis of IBD and IBS in complex trait studies. Nature Reviews Genetics, 11(11), 800-805.}
#'
#' @export

Gmatrix <- function (SNPmatrix = NULL, method = 2, missingValue = -9, maf = 0, verify.posdef = TRUE) 
{
    Time = proc.time()
    if (is.null(SNPmatrix)) {
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
    if(class(SNPmatrix)!="matrix"){
        cat("SNPmatrix class is:",class(SNPmatrix),"\n")
        stop("SNPmatrix class must be matrix. Please verify it.")
        }
    if (!is.na(missingValue)) {
        m <- match(SNPmatrix, missingValue, 0)
        SNPmatrix[m > 0] <- NA
    }
    NumberMarkers = ncol(SNPmatrix)
    nindTotal = colSums(!is.na(SNPmatrix))
    nindAbs = max(nindTotal)
    cat("\tNumber of Individuals:", max(nindTotal), "\n")
    cat("\tNumber of Markers:", NumberMarkers, "\n")
    alelleFreq = function(x, y) {
        (2 * length(which(x == y)) + length(which(x == 1)))/(2 * 
            length(which(!is.na(x))))
    }
    Frequency = cbind(apply(SNPmatrix, 2, function(x) alelleFreq(x, 
        0)), apply(SNPmatrix, 2, function(x) alelleFreq(x, 2)))
    if (any(Frequency[, 1] <= maf) & maf != 0) {
        cat("\t", length(which(Frequency[, 1] <= maf)), "markers dropped due to maf cutoff of", maf, "\n")
        SNPmatrix = SNPmatrix[,-which(Frequency[, 1] <= maf)]
        cat("\t", ncol(SNPmatrix), "markers kept \n")
        Frequency = as.matrix(Frequency[-which(Frequency[,1] <= 
            maf), ])
        NumberMarkers = ncol(SNPmatrix)
    }
    Sum2pq = 2 * t(Frequency[, 1]) %*% Frequency[, 2]
    FreqP = matrix(rep(Frequency[, 2], each = nrow(SNPmatrix)), 
        ncol = ncol(SNPmatrix))
    if (method == "MarkersMatrix") {
        Gmatrix <- !is.na(SNPmatrix)
        Gmatrix <- tcrossprod(Gmatrix, Gmatrix)
        return(Gmatrix)
    }
    if (method == "VanRaden") {
       # SNP.VanRaden = SNPmatrix
        SNPmatrix = SNPmatrix- 2 * FreqP
        SNPmatrix[is.na(SNPmatrix)] <- 0
        Gmatrix = (tcrossprod(SNPmatrix, SNPmatrix))/as.numeric(Sum2pq)
    }
    if (method == "Powell") {
        FreqPQ = matrix(rep(2 * Frequency[, 1] * Frequency[, 
            2], each = nrow(SNPmatrix)), ncol = ncol(SNPmatrix))
        G.all = (SNPmatrix^2 - (1 + 2 * FreqP) * SNPmatrix + 
            2 * FreqP^2)/FreqPQ
        G.ii = as.matrix(colSums(t(G.all), na.rm = T))
        SNPmatrix = (SNPmatrix - (2 * FreqP))/sqrt(FreqPQ)
        G.ii.hat = 1 + (G.ii)/NumberMarkers
        SNPmatrix[is.na(SNPmatrix)] <- 0
        Gmatrix = (tcrossprod(SNPmatrix, SNPmatrix))/NumberMarkers
        diag(Gmatrix) = G.ii.hat
    }
    if (verify.posdef) {
        e.values <- eigen(Gmatrix, symmetric = TRUE)$values
        indicator <- length(which(e.values <= 0))
        if (indicator > 0) 
            cat("\t Matrix is NOT posiive definite. It has ", indicator, 
                " eigenvalues <= 0 \n \n")
    }
    Time = as.matrix(proc.time() - Time)
    cat("Completed! Time =", Time[3], " seconds \n")
    return(Gmatrix)
}

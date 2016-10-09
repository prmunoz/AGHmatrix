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
# Last update: 08-Oct-2016 						#
# License: GNU General Public License version 2 (June, 1991) or later 	#
# 									#
#########################################################################

#' Construction of Relationship Matrix G
#'
#' Given a matrix (individual x markers), a method, a missing value, and a maf threshold, return a additive or dominance relationship matrix.
#'
#' @param SNPmatrix matrix (n x m), where n is is individual names and m is marker names (coded inside the matrix as 0, 1, 2, missingValue). 
#' @param method "Yang" or "VanRaden" for marker-based additive relationship matrix. "Su" or "Vitezica" for marker-based dominance relationship matrix. "MarkersMatrix" for a matrix with the amount of shared markers between individuals (3). Default="Yang".
#' @param missingValue missing value in data. Default=-9.
#' @param maf max of missing data accepted to each marker. Default=0.05.
#' @param verify.posdef verify if the resulting matrix is positive-definite. Default=TRUE.
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
#' @author Rodrigo R Amadeu \email{rramadeu@@gmail.com} and Marcio Resende Jr.
#' 
#' @references \emph{Su, G., Christensen, O.F., Ostersen, T., Henryon, M. and Lund, M.S., 2012. Estimating additive and non-additive genetic variances and predicting genetic merits using genome-wide dense single nucleotide polymorphism markers. PloS one, 7(9), p.e45293.}
#' @references \emph{VanRaden, P.M., 2008. Efficient methods to compute genomic predictions. Journal of dairy science, 91(11), pp.4414-4423.}
#' @references \emph{Vitezica, Z.G., Varona, L. and Legarra, A., 2013. On the additive and dominant variance and covariance of individuals within the genomic selection scope. Genetics, 195(4), pp.1223-1230.}
#' @references \emph{Yang, J., Benyamin, B., McEvoy, B.P., Gordon, S., Henders, A.K., Nyholt, D.R., Madden, P.A., Heath, A.C., Martin, N.G., Montgomery, G.W. and Goddard, M.E., 2010. Common SNPs explain a large proportion of the heritability for human height. Nature genetics, 42(7), pp.565-569.}
#'
#' @export

Gmatrix <- function (SNPmatrix = NULL, method = "Yang", missingValue = -9, maf = 0, verify.posdef = TRUE) 
{
    Time = proc.time()
    if (is.null(SNPmatrix)) {
        stop(deparse("Please define the variable SNPdata"))
    }
    if (all(method != c("Yang", "VanRaden", "Su", "Vitezica", "MarkersMatrix"))) {
        stop("Method to build Gmatrix has to be either `Yang` or `VanRaden` for marker-based additive relationship matrix, or
`Su` or `Vitezica` for marker-based dominance relationship matrx, or
 `MarkersMatrix` for the a matrix with amount of shared-marks by individuals pairs")
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
    FreqP = matrix(rep(Frequency[, 2], each = nrow(SNPmatrix)), 
        ncol = ncol(SNPmatrix))
    if (method == "MarkersMatrix") {
        Gmatrix <- !is.na(SNPmatrix)
        Gmatrix <- tcrossprod(Gmatrix, Gmatrix)
        return(Gmatrix)
    }
    if (method == "VanRaden") {
        Sum2pq = 2 * crossprod(Frequency[,1],Frequency[,2])
        SNPmatrix = SNPmatrix- 2 * FreqP
        SNPmatrix[is.na(SNPmatrix)] <- 0
        Gmatrix = (tcrossprod(SNPmatrix, SNPmatrix))/as.numeric(Sum2pq)
    }
    if (method == "Yang") {
        FreqPQ = matrix(rep(2 * Frequency[, 1] * Frequency[,2], each = nrow(SNPmatrix)), ncol = ncol(SNPmatrix))
        G.all = (SNPmatrix^2 - (1 + 2 * FreqP) * SNPmatrix + 
            2 * FreqP^2)/FreqPQ
        G.ii = as.matrix(colSums(t(G.all), na.rm = T))
        SNPmatrix = (SNPmatrix - (2 * FreqP))/sqrt(FreqPQ)
        G.ii.hat = 1 + (G.ii)/NumberMarkers
        SNPmatrix[is.na(SNPmatrix)] <- 0
        Gmatrix = (tcrossprod(SNPmatrix, SNPmatrix))/NumberMarkers
        diag(Gmatrix) = G.ii.hat
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
        SNPmatrix[is.na(SNPmatrix)] <- -9
        SNPmatrix <- (SNPmatrix==0)*-2*(FreqP^2) +
            (SNPmatrix==1)*2*(FreqP)*(1-FreqP) +
            (SNPmatrix==2)*-2*((1-FreqP)^2)
        Gmatrix <- tcrossprod(SNPmatrix,SNPmatrix)/sum(TwoPQ^2)
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
    return(Gmatrix)
}

#####################################################################
# 									
# Package: AGHmatrix 							
# 									
# File: Gmatrix.R
# Contains: Gmatrix slater_par check_Gmatrix_data			
# 									
# Written by Rodrigo Rampazo Amadeu 			
# Contributors: Marcio Resende Jr, Leticia AC Lara, Ivone Oliveira, Luis Felipe V Ferrao
# 									
# First version: Feb-2014 					
# Last update: 05-Aug-2021 						
# License: GPL-3	
# 									
#####################################################################

#' Construction of Relationship Matrix G
#'
#' Given a matrix (individual x markers), a method, a missing value, and a maf 
#' threshold, return a additive or non-additive relationship matrix.
#' For diploids, the methods "Yang" and "VanRaden" for additive relationship 
#' matrices, and "Su" and "Vitezica" for non-additive relationship matrices 
#' are implemented. For autopolyploids, the method "VanRaden" for additive 
#' relationship, method "Slater" for full-autopolyploid model including 
#' non-additive effects, and pseudo-diploid parametrization are implemented.
#' Weights are implemented for "VanRaden" method as described in Liu (2020).
#' 
#' @param SNPmatrix matrix (n x m), where n is is individual names and m is 
#' marker names (coded inside the matrix as 0, 1, 2, ..., ploidy, and, 
#' missingValue). 
#' @param method "Yang" or "VanRaden" for marker-based additive relationship 
#' matrix. "Su" or "Vitezica" for marker-based dominance relationship matrix. 
#' "Slater" for full-autopolyploid model including non-additive effects.
#' "Endelman" for autotetraploid dominant (digentic) relationship matrix. 
#' "MarkersMatrix" for a matrix with the amount of shared markers between 
#' individuals (3).
#' Default is "VanRaden", for autopolyploids will be computed a scaled 
#' product (similar to Covarrubias-Pazaran, 2006).
#' @param missingValue missing value in data. Default=-9.
#' @param thresh.missing threshold on missing data, SNPs below of this 
#' frequency value will be maintained, if equal to 1, no threshold and 
#' imputation is considered. Default = 1.
#' @param maf minimum allele frequency accepted to each marker. Default=0.
#' @param verify.posdef verify if the resulting matrix is positive-definite. 
#' Default=FALSE.
#' @param ploidy data ploidy (an even number between 2 and 20). Default=2.
#' @param pseudo.diploid if TRUE, uses pseudodiploid parametrization of 
#' Slater (2016).
#' @param ratio if TRUE, molecular data are considered ratios and its 
#' computed the scaled product of the matrix (as in "VanRaden" method).
#' @param impute.method "mean" to impute the missing data by the mean per 
#' marker, "mode" to impute the missing data by the mode per marker, 
#' "global.mean" to impute the missing data by the mean across all markers,
#' "global.mode" to impute the missing data my the mode across all marker. 
#' Default = "none".
#' @param integer if FALSE, not check for integer numbers. Default=TRUE.
#' @param ratio.check if TRUE, run Mcheck with ratio data.
#' @param weights vector with weights for each marker. Only works if 
#' method="VanRaden". Default is a vector of 1's (equal weight).
#' @param ploidy.correction It sets the denominator (correction) of the 
#' rossprod. Used only when ploidy > 2 for "VanRaden" and ratio models.
#' If TRUE, it uses the sum of "Ploidy" times "Frequency" times "(1-Frequency)" 
#' of each marker as method 1 in VanRaden 2008 and Endelman (2018).
#' When ratio=TRUE, it uses "1/Ploidy" times "Frequency" times "(1-Frequency)". 
#' If FALSE, it uses the sum of the sampling variance of each marker. 
#' Default = FALSE. 
#' @param rmv.mono if monomorphic markers should be removed. Default=FALSE.
#' @param thresh.htzy threshold heterozigosity, remove SNPs below this 
#' threshold. Default=0.
#' @param ASV if TRUE, transform matrix into average semivariance (ASV) 
#' equivalent (K = K / (trace(K) / (nrow(K)-1))). Details formula 2 of 
#' Fieldmann et al. (2022). Default = FALSE.
#' @return Matrix with the marker-bases relationships between the individuals
#'
#' @useDynLib AGHmatrix, .registration = TRUE
#' @importFrom Rcpp evalCpp
#' @importFrom matrixStats colVars
#' @examples
#' \dontrun{
#' ## Diploid Example
#' data(snp.pine)
#' str(snp.pine)
#' Gmatrix.Yang <- Gmatrix(snp.pine, method="Yang", missingValue=-9, maf=0.05)
#' Gmatrix.VanRaden <- Gmatrix(snp.pine, method="VanRaden", missingValue=-9, maf=0.05)
#' Gmatrix.Su <- Gmatrix(snp.pine, method="Su", missingValue=-9, maf=0.05)
#' Gmatrix.Vitezica <- Gmatrix(snp.pine, method="Vitezica", missingValue=-9, maf=0.05)
#' 
#' ## Autetraploid example
#' data(snp.sol)
#' Gmatrix.VanRaden <- Gmatrix(snp.sol, method="VanRaden", ploidy=4)
#' Gmatrix.Endelman <- Gmatrix(snp.sol, method="Endelman", ploidy=4) 
#' Gmatrix.Slater <- Gmatrix(snp.sol, method="Slater", ploidy=4)
#' Gmatrix.Pseudodiploid <- Gmatrix(snp.sol, method="VanRaden", ploidy=4, pseudo.diploid=TRUE) 
#' Gmatrix.weighted <- Gmatrix(snp.sol, method="VanRaden", weights = runif(3895,0.001,0.1), ploidy=4)
#' }
#' 
#' @author Rodrigo R. Amadeu \email{rramadeu@gmail.com}
#' @author Thiago de Paula Oliveira \email{toliveira@abacusbio.com}
#' @author Marcio Resende Jr.
#' @author Letícia A. C. Lara
#' @author Ivone Oliveira
#' @author Felipe V. Ferrao
#' @references \emph{Covarrubias-Pazaran, G. 2016. Genome assisted prediction 
#' of quantitative traits using the R package sommer. PLoS ONE 11(6):1-15.}
#' @references \emph{Endelman, JB, et al., 2018. Genetic variance partitioning 
#' and genome-wide prediction with allele dosage information in autotetraploid 
#' potato. Genetics, 209(1) pp. 77-87.}
#' @references \emph{Feldmann MJ, et al. 2022. Average semivariance directly 
#' yields accurate estimates of the genomic variance in complex trait analyses. 
#' G3 (Bethesda), 12(6).}
#' @references \emph{Liu, A, et al. 2020. Weighted single-step genomic best 
#' linear unbiased prediction integrating variants selected from sequencing 
#' data by association and bioinformatics analyses. Genet Sel Evol 52, 48.}
#' @references \emph{Slater, AT, et al. 2016. Improving genetic gain with 
#' genomic selection in autotetraploid potato. The Plant Genome 9(3), pp.1-15.}
#' @references \emph{Su, G, et al. 2012. Estimating additive and non-additive 
#' genetic variances and predicting genetic merits using genome-wide dense 
#' single nucleotide polymorphism markers. PloS one, 7(9), p.e45293.}
#' @references \emph{VanRaden, PM, 2008. Efficient methods to compute genomic 
#' predictions. Journal of dairy science, 91(11), pp.4414-4423.}
#' @references \emph{Vitezica, ZG, et al. 2013. On the additive and dominant 
#' variance and covariance of individuals within the genomic selection scope. 
#' Genetics, 195(4), pp.1223-1230.}
#' @references \emph{Yang, J, et al. 2010. Common SNPs explain a large 
#' proportion of the heritability for human height. Nature genetics, 42(7), 
#' pp.565-569.}
#'
#' @export
Gmatrix <- function(SNPmatrix = NULL, 
                    method = "VanRaden", 
                    missingValue = -9, 
                    maf = 0, 
                    thresh.missing = 1,
                    verify.posdef = FALSE, 
                    ploidy = 2,
                    pseudo.diploid = FALSE, 
                    integer = TRUE,
                    ratio = FALSE, 
                    impute.method = "none", 
                    rmv.mono = FALSE, 
                    thresh.htzy = 0,
                    ratio.check = TRUE,
                    weights = NULL, 
                    ploidy.correction = FALSE, 
                    ASV = FALSE) {
  Time <- proc.time()
  #-----------------------------------------------------------------------------
  # ==== Initial preprocessing ====
  #-----------------------------------------------------------------------------
  check_matrix((SNPmatrix))
  SNPmatrix <- as.matrix(data.matrix(SNPmatrix))
  ids <- rownames(SNPmatrix)
  
  if (ploidy %% 2 != 0 || ploidy < 2 || ploidy > 20) {
    stop("ploidy must be an even number between 2 and 20")
  }
  
  markers <- colnames(SNPmatrix)
  
  if (!is.null(weights) && length(weights) != ncol(SNPmatrix)) {
    stop("weights length (", length(weights), 
         ") does not match number of markers (", ncol(SNPmatrix), ")")
  }
  
  if (ratio) method <- "VanRaden"
  
  if (!is.na(missingValue)) {
    SNPmatrix[SNPmatrix == missingValue] <- NA
  }
  
  check_Gmatrix_data(SNPmatrix = SNPmatrix, 
                     method = method, 
                     ploidy = ploidy, 
                     ratio = ratio, 
                     integer = integer)
  
  # Early return for MarkersMatrix using the raw NA pattern (parity with legacy)
  if (identical(method, "MarkersMatrix")) {
    if (exists("Gmatrix_MarkersMask", mode = "function")) {
      Gmatrix <- Gmatrix_MarkersMask(SNPmatrix)
    } else {
      mask <- !is.na(SNPmatrix)
      storage.mode(mask) <- "double"
      Gmatrix <- tcrossprod(mask, mask)
    }
    if (!is.null(ids)) dimnames(Gmatrix) <- list(ids, ids)
    return(Gmatrix)
  }
  
  #-----------------------------------------------------------------------------
  # ==== Missing-data checks / filtering (only when requested) ====
  #-----------------------------------------------------------------------------
  do_mcheck <- (!ratio || (ratio && ratio.check)) &&
    !(isTRUE(all.equal(thresh.missing, 1)) || identical(impute.method, "none"))
  
  if (do_mcheck) {
    SNPmatrix <- Mcheck(
      SNPmatrix,
      ploidy          = ploidy,
      thresh.maf      = maf,
      rmv.mono        = rmv.mono,
      thresh.htzy     = thresh.htzy,
      thresh.missing  = thresh.missing,
      impute.method   = impute.method
    )
  }
  
  #-----------------------------------------------------------------------------
  # ==== Diagnostics ====
  #-----------------------------------------------------------------------------
  NumberMarkers <- ncol(SNPmatrix)
  nindTotal <- colSums(!is.na(SNPmatrix))
  message("Initial data: ")
  message("\tNumber of Individuals:", max(nindTotal))
  message("\tNumber of Markers:", NumberMarkers)
  message("Building G matrix using method: ", method, ", ploidy: ", ploidy)
  
  #-----------------------------------------------------------------------------
  ## ===== Prepare frequencies =====
  #-----------------------------------------------------------------------------
  P <- colMeans(SNPmatrix, na.rm = TRUE)
  if (ploidy == 2) { 
    fr   <- diploid_p0p2_TwoPQ_cpp(SNPmatrix)
    p0   <- as.numeric(fr$p0)
    p2   <- as.numeric(fr$p2)
    TwoPQ <- as.numeric(fr$TwoPQ)
    
    Frequency <- cbind(p0, p2)
    FreqP     <- matrix(p2, nrow = nrow(SNPmatrix), ncol = ncol(SNPmatrix), 
                        byrow = TRUE)
    }
  
  if (ploidy > 2 && pseudo.diploid) {
    SNPmatrix[, P > ploidy / 2] <- ploidy - SNPmatrix[, P > ploidy / 2]
    Frequency <- colMeans(SNPmatrix, na.rm = TRUE) / ploidy
    Frequency <- cbind(1 - Frequency, Frequency)
    FreqP <- matrix(Frequency[, 2], nrow = nrow(SNPmatrix), 
                    ncol = ncol(SNPmatrix), byrow = TRUE)
    SNPmatrix[SNPmatrix %in% c(1:(ploidy - 1))] <- 1
    SNPmatrix[SNPmatrix == ploidy] <- 2
  }
  
  #-----------------------------------------------------------------------------
  ## ===== Main method  =====
  #-----------------------------------------------------------------------------
  if (method == "VanRaden") {
    if (is.null(weights)) {
      if (ploidy == 2 && !ratio) {
        Gmatrix <- Gmatrix_vanraden(SNPmatrix, Frequency[, 2], TwoPQ)
      } else {
        Gmatrix <- Gmatrix_vanraden_poly_unweighted(SNPmatrix, ploidy, ratio, 
                                                    ploidy.correction)
      }
    } else {
      weights <- weights[match(colnames(SNPmatrix), markers)]
      if (ploidy == 2 && !ratio) {
        Gmatrix <- Gmatrix_VanRaden_weighted(SNPmatrix, weights, 
                                             Frequency[, 2], TwoPQ)
      } else {
        Gmatrix <- 
          Gmatrix_vanraden_poly_weighted(SNPmatrix, as.numeric(weights), 
                                         ploidy, ratio, ploidy.correction)
      }
    }
  }
  
  if (method == "Yang") {
    if (ploidy != 2) stop("Yang method is defined for diploids (ploidy = 2).")
    FreqPQ <- matrix(rep(2 * Frequency[, 1] * Frequency[, 2], 
                         each = nrow(SNPmatrix)),
                     ncol = ncol(SNPmatrix))
    G.all <- (SNPmatrix^2 - (1 + 2 * FreqP) * SNPmatrix + 2 * (FreqP^2)) / FreqPQ
    G.ii <- as.matrix(colSums(t(G.all), na.rm = TRUE))
    Z <- (SNPmatrix - (2 * FreqP)) / sqrt(FreqPQ)
    G.ii.hat <- 1 + (G.ii) / NumberMarkers
    Z[is.na(Z)] <- 0
    Gmatrix <- (tcrossprod(Z, Z)) / NumberMarkers
    diag(Gmatrix) <- G.ii.hat
  }
  
  if (method == "Su") {
    Gmatrix <- Gmatrix_Su(SNPmatrix)
  }
  
  if (method == "Vitezica") {
    Gmatrix <- Gmatrix_Vitezica(SNPmatrix, FreqP)
  }
  
  if (method == "Slater") {
    P <- colMeans(SNPmatrix, na.rm = TRUE)
    SNPmatrix[, which(P > ploidy/2)] <- ploidy - SNPmatrix[, which(P > ploidy/2)]
    SNPmatrix <- slater_par_cpp(SNPmatrix, ploidy = ploidy)
    
    NumberMarkers <- ncol(SNPmatrix)
    Frequency     <- colMeans(SNPmatrix, na.rm = TRUE)
    FreqP         <- matrix(rep(Frequency, each = nrow(SNPmatrix)),
                            nrow = nrow(SNPmatrix))
    
    drop.alleles  <- which(Frequency == 0)
    if (length(drop.alleles) > 0) {
      Frequency <- Frequency[-drop.alleles]
      SNPmatrix <- SNPmatrix[, -drop.alleles, drop = FALSE]
      FreqP     <- FreqP[, -drop.alleles, drop = FALSE]
    }
    
    Gmatrix <- Gmatrix_Slater(SNPmatrix, Frequency, FreqP, NumberMarkers) 
  }
  
  if (method == "Endelman") {
    if (ploidy != 4) stop("'Endelman' method is just implemented for ploidy=4")
    Gmatrix <- Gmatrix_Endelman(SNPmatrix, ploidy)
  }
  
  if (verify.posdef) {
    e.values <- eigen(Gmatrix, symmetric = TRUE)$values
    indicator <- sum(e.values <= 0)
    if (indicator > 0) message("\t Matrix is NOT positive definite. It has ", 
                               indicator, " eigenvalues <= 0")
  }
  
  if (ASV) {
    Gmatrix <- get_ASV(Gmatrix)
  }
  if (!is.null(ids)) dimnames(Gmatrix) <- list(ids, ids)
  
  message("Completed! Time = ", round(proc.time()[3] - Time[3], 2), " seconds")
  
  attr(Gmatrix, "method") <- method
  attr(Gmatrix, "ploidy") <- ploidy
  attr(Gmatrix, "nmarkers") <- ncol(SNPmatrix)
  return(Gmatrix)
}


## Internal Functions ##
get_ASV = function(x){
  return( x / ( sum(diag(x)) / (nrow(x) - 1)) )
}

# Internal function to check input Gmatrix arguments
check_Gmatrix_data <- function(SNPmatrix,ploidy,method, ratio=FALSE, integer=TRUE){
  if (is.null(SNPmatrix)) {
    stop(deparse("Please define the variable SNPdata"))
  }
  if (all(method != c("Yang", "VanRaden", "Slater", "Su", "Vitezica", "MarkersMatrix","Endelman"))) {
    stop("Method to build Gmatrix has to be either `Yang` or `VanRaden` for marker-based additive relationship matrix, or `Su` or `Vitezica` or `Endelman` for marker-based dominance relationship matrx, or `MarkersMatrix` for matrix with amount of shared-marks by individuals pairs")
  }
  
  #  if( method=="Yang" && ploidy>2)
  #    stop("Change method to 'VanRaden' for ploidies higher than 2 for marker-based additive relationship matrix")
  
  if( method=="Su" && ploidy>2)
    stop("Change method to 'Slater' for ploidies higher than 2 for marker-based non-additive relationship matrix")
  
  if( method=="Vitezica" && ploidy>2)
    stop("Change method to 'Slater' for ploidies higher than 2 for marker-based non-additive relationship matrix")
  
  if(!is.matrix(SNPmatrix)){
    cat("SNPmatrix class is:",class(SNPmatrix),"\n")
    stop("SNPmatrix class must be matrix. Please verify it.")
  }
  
  if(!ratio){
    if( ploidy > 20 | (ploidy %% 2) != 0)
      stop(deparse("Only even ploidy from 2 to 20"))
    
    t <- max(SNPmatrix,na.rm = TRUE)
    if( t > ploidy )
      stop(deparse("Check your data, it has values above ploidy number"))
    
    t <- min(SNPmatrix,na.rm=TRUE)
    if( t < 0 )
      stop(deparse("Check your data, it has values under 0"))
    
    if(integer)
      if(prod(SNPmatrix == round(SNPmatrix),na.rm = TRUE)==0)
        stop(deparse("Check your data, it has not integer values"))
  }
  
  if(ratio){
    t <- max(SNPmatrix,na.rm = TRUE)
    if( t > 1)
      stop(deparse("Check your data, it has values above 1. It is expected a ratio values [0;1]."))
    
    t <- min(SNPmatrix,na.rm=TRUE)
    if( t < 0 )
      stop(deparse("Check your data, it has values under 0. It is expected a ratio values [0;1]."))
  }
}

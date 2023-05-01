#########################################
# 									
# Package: AGHmatrix 							
# 									
# File: Mcheck.R							
# Contains: Mcheck							
# 									
# Written by Luis F V Ferrao and Rodrigo Amadeu		
# 									
# First version: Feb-2014 						
# Last update: 30-Mar-2023 						
# License: GPL-3
# 									
#########################################

#' Check and filter markers
#'
#' This function does different filtering on the marker matrix
#'
#' @param SNPmatrix matrix (n x m), where n is is individual names and m is marker names (coded inside the matrix as 0, 1, 2, ..., ploidy, and, missingValue). 
#' @param ploidy data ploidy (an even number between 2 and 20). Default=2.
#' @param missingValue missing value in data. Default=-9.
#' @param thresh.missing threshold on missing data, SNPs below of this frequency value will be maintained, if equal to 1, no threshold and imputation is considered. Default = 0.50.
#' @param thresh.maf minimum allele frequency accepted to each marker. Default=0.05.
#' @param thresh.htzy threshold heterozigosity, remove SNPs below this threshold. Default=0.
#' @param impute.method "mean" to impute the missing data by the mean per marker, "mode" to impute the missing data by the mode per marker, "global.mean" to impute the missing data by the mean across all markers, "global.mode" to impute the missing data my the mode across all marker. Default = "mean".
#' @param rmv.mono if monomorphic markers should be removed. Default=TRUE.
#'
#' @return SNPmatrix after filtering steps.
#'
#' @examples 
#' data(snp.pine)
#' M = Mcheck(snp.pine)
#'
#' @author Luis F V Ferrao and Rodrigo Amadeu, \email{rramadeu@@gmail.com}
#'
#' @export
#' 
Mcheck = function(SNPmatrix = NULL,
                  ploidy=2,
                  missingValue = -9,
                  thresh.maf = 0.05,
                  thresh.missing = 0.9,
                  thresh.htzy = 0, 
                  impute.method = "mean",
                  rmv.mono=TRUE){
  
  # SNP missing data
  ncol.init <- ncol(SNPmatrix)
  
  if (!is.na(missingValue)) {
    m <- match(SNPmatrix, missingValue, 0)
    SNPmatrix[m > 0] <- NA
  }
  
  missing <- apply(SNPmatrix, 2, function(x) sum(is.na(x))/nrow(SNPmatrix))
  missing.low = missing <= thresh.missing
  cat("\nMissing data check: \n")
  if(any(missing.low)){
    cat("\tTotal SNPs:", ncol(SNPmatrix),"\n")
    cat("\t",ncol(SNPmatrix) - sum(missing.low), "SNPs dropped due to missing data threshold of", thresh.missing,"\n")
    cat("\tTotal of:",sum(missing.low), " SNPs \n")
    idx.rm <- which(missing.low)
    SNPmatrix <- SNPmatrix[, idx.rm, drop=FALSE]
  } else{
    cat("\tNo SNPs with missing data, missing threshold of = ", thresh.missing,"\n")
  }
  
  # Minor alele frequency
  MAF <- apply(SNPmatrix, 2, function(x) {
    AF <- mean(x, na.rm = T)/ploidy
    MAF <- ifelse(AF > 0.5, 1 - AF, AF) # Minor allele freq can be ref allele or not
  })
  snps.low <- MAF < thresh.maf
  cat("\nMAF check: \n")
  if(any(snps.low)){
    cat("\t",sum(snps.low), "SNPs dropped with MAF below", thresh.maf,"\n")
    cat("\tTotal:",ncol(SNPmatrix) - sum(snps.low), " SNPs \n")
    idx.rm <- which(snps.low)
    SNPmatrix <- SNPmatrix[, -idx.rm, drop=FALSE]
  } else{
    cat("\tNo SNPs with MAF below", thresh.maf,"\n")
  }
  
  # monomorphic SNPs
  if(rmv.mono){
    mono <- apply(SNPmatrix, 2, function(x) {
      equal <- isTRUE(all.equal(x, rep(x[1], length(x))))
    })
    cat("\nMonomorphic check: \n")
    if(any(mono)){
      cat("\t",sum(mono), "monomorphic SNPs \n")
      cat("\tTotal:",ncol(SNPmatrix) - sum(mono), "SNPs \n")
      idx.rm <- which(mono)
      SNPmatrix <- SNPmatrix[, -idx.rm, drop=FALSE]
    } else{
      cat("\tNo monomorphic SNPs \n")
    }
  }
  
  
  # Imputation
  if(impute.method=="global.mean"){
    ix <- which(is.na(SNPmatrix))
    if (length(ix) > 0) {
      SNPmatrix[ix] <- mean(SNPmatrix,na.rm = TRUE)
    }
  }
  
  if(impute.method=="global.mode"){
    ix <- which(is.na(SNPmatrix))
    if (length(ix) > 0) {
      SNPmatrix[ix] <- as.integer(names(which.max(table(SNPmatrix))))
    }
  }
  
  if(impute.method=="mean"){
    imputvalue = apply(SNPmatrix,2,mean,na.rm=TRUE)
    ix = which(is.na(SNPmatrix),arr.ind=TRUE)
    SNPmatrix[ix] = imputvalue[ix[,2]]
  }
  
  if(impute.method=="mode"){
    imputvalue = apply(SNPmatrix, 2, function(x) as.integer(names(which.max(table(x)))))
    ix = which(is.na(SNPmatrix),arr.ind=TRUE)
    SNPmatrix[ix] = imputvalue[ix[,2]]
  }    
  
  
  # Heterozigosity
  htrz <- apply(SNPmatrix, 2, function(x) sum( x!= 0 & x != ploidy,na.rm=T)/nrow(SNPmatrix))
  htrz.low = htrz < thresh.htzy 
  cat("\nHeterozigosity data check: \n")
  if(any(htrz.low)){
    cat("\tTotal SNPs:", ncol(SNPmatrix),"\n")
    cat("\t",ncol(SNPmatrix) - sum(htrz.low), "SNPs dropped due to heterozygosity threshold of", thresh.htzy,"\n")
    cat("\tTotal of:",sum(htrz.low), " SNPs \n")
    idx.rm <- which(htrz.low)
    SNPmatrix <- SNPmatrix[, idx.rm, drop=FALSE]
  } else{
    cat("\tNo SNPs with heterozygosity, missing threshold of = ", thresh.htzy,"\n")
  }
  
  # Total of SNPs
  cat("\nSummary check: \n")
  cat("\tInitial: ", ncol.init, "SNPs \n")
  cat("\tFinal: ", ncol(SNPmatrix), " SNPs (", ncol.init - ncol(SNPmatrix), " SNPs removed) \n \n")
  return(SNPmatrix)
}

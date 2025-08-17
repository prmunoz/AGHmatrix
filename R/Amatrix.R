#########################################################################
# 									
# Package: AGHmatrix 							
# 									
# File: Amatrix.R
# Contains: Amatrix 							
# 									
# Written by Rodrigo Rampazo Amadeu and Thiago de Paula Oliveira
# 									
# First version: Feb-2014 
# Last update: 12-July-2024 						
# License: GPL-3
# 								
#########################################################################

#' @title Construction of Relationship Matrix A
#'
#' @description
#' Constructs an additive relationship matrix (A) from pedigree data in a 
#' three-column format. The calculation depends on the specified ploidy level 
#' (must be an even  number). If \code{ploidy = 4}, the matrix can incorporate a 
#' user-defined proportion (\code{w}) of parental gametes that are IBD 
#' (Identical by Descent) due to double reduction.
#'
#' When \code{ploidy = 2} and \code{dominance = FALSE}, the diploid additive 
#' numerator relationship matrix is built using Henderson's (1976) method. If 
#' \code{dominance = TRUE}, a diploid dominance numerator relationship matrix 
#' is computed following Cockerham (1954). The recursive approach is detailed 
#' in Mrode (2005).
#'
#' For \code{ploidy > 2} and \code{dominance = FALSE}, an autopolyploid 
#' additive relationship matrix is computed based on Kerr et al. (2012). 
#' If \code{slater = TRUE} and \code{ploidy = 4}, Slater's (2013) method is 
#' used instead.
#'
#' @param data Pedigree data frame in a 3-column format. Unknown values should 
#' be set to 0.
#' @param ploidy An even number representing the ploidy level (default is 2).
#' @param w Proportion of parental gametes that are IBD due to double 
#' reduction (default is 0). Applicable only when \code{ploidy = 4}.
#' @param verify Logical; if TRUE (default), the pedigree file is checked for 
#' conflicting entries.
#' @param dominance Logical; if TRUE, computes the dominance relationship 
#' matrix (valid only for \code{ploidy = 2}).
#' @param slater Logical; if TRUE and \code{ploidy = 4}, computes the additive 
#' autotetraploid relationship matrix using Slater (2013).
#' @param ASV Logical; if TRUE, transforms the matrix into the average 
#' semivariance (ASV) defined as 
#' \eqn{K = K / (\text{trace}(K) / (\text{nrow}(K) - 1))}. See Feldmann et al. 
#' (2022), Formula 2.
#' @param ... Additional arguments passed to \code{datatreat()}.
#' 
#' @return A square matrix representing the relationship between individuals 
#' in the pedigree.
#'
#' @examples
#' data(ped.mrode)
#' # Additive relationship matrix for diploids (Henderson 1976):
#' Amatrix(ped.mrode, ploidy = 2)
#'
#' # Dominance relationship matrix for diploids (Cockerham 1954):
#' Amatrix(ped.mrode, ploidy = 2, dominance = TRUE)
#'
#' # Additive relationship matrix for autotetraploids (Kerr 2012):
#' Amatrix(ped.mrode, ploidy = 4)
#'
#' # Additive relationship matrix for autooctaploids (Kerr 2012):
#' Amatrix(ped.mrode, ploidy = 8)
#'
#' # Additive relationship matrix for autotetraploids with double reduction 
#' # (Kerr 2012):
#' Amatrix(ped.mrode, ploidy = 4, w = 0.1)
#'
#' # Additive relationship matrix for autotetraploids with double reduction 
#' # using Slater (2013):
#' Amatrix(ped.mrode, ploidy = 4, w = 0.1, slater = TRUE)
#'
#' # Additive relationship matrix for autohexaploids with double reduction 
#' # (Kerr 2012):
#' Amatrix(ped.mrode, ploidy = 6, w = 0.1)
#'
#' @useDynLib AGHmatrix, .registration = TRUE
#' @importFrom Rcpp evalCpp
#' @note This function uses compiled C++ code via Rcpp for improved performance.
#'
#' @author Rodrigo R Amadeu, \email{rramadeu@@gmail.com}
#' @author Thiago de Paula Oliveira, \email{toliveira@@abacusbio.com}
#'
#' @references 
#'   \emph{Cockerham, C.C. (1954). An extension of the concept of partitioning 
#'   hereditary variance for analysis of covariances among relatives when 
#'   epistasis is present. Genetics, 39, 859–882.}
#' @references 
#'   \emph{Feldmann, M.J., et al. (2022). Average semivariance directly yields 
#'   accurate estimates of the genomic variance in complex trait analyses. 
#'   G3 (Bethesda), 12(6).}
#' @references 
#'   \emph{Henderson, C.R. (1976). A simple method for computing the inverse 
#'   of a numerator relationship matrix used in prediction of breeding values. 
#'   Biometrics, 32, 69–83.}
#' @references 
#'   \emph{Kerr, R.J., et al. (2012). Use of the numerator relationship matrix 
#'   in genetic analysis of autopolyploid species. Theoretical and Applied 
#'   Genetics, 124, 1271–1282.}
#' @references 
#'   \emph{Mrode, R.A. (2014). Linear Models for the Prediction of Animal 
#'   Breeding Values (3rd ed.). CABI.}
#' @references 
#'   \emph{Slater, A.T., et al. (2013). Improving the analysis of low 
#'   heritability complex traits for enhanced genetic gain in potato. 
#'   Theoretical and Applied Genetics, 127, 809–820.}
Amatrix <- function(data = NULL,
                    ploidy = 2,
                    w = 0,
                    verify = TRUE,
                    dominance = FALSE,
                    slater = FALSE,
                    ASV = FALSE,
                    ...) {
  #-------------------------------------------------------------------------
  # 1. Validate Inputs
  #-------------------------------------------------------------------------
  if (is.null(data) || !is.data.frame(data) || ncol(data) < 3)
    stop("Input 'data' must be a non-null data frame with at least 3 columns.")
  
  if (nrow(data) == 0)
    stop("Pedigree data is empty.")
  
  if (ploidy %% 2 != 0)
    stop("Ploidy should be an even number.")
  
  if (ploidy < 2)
    stop("Ploidy must be at least 2.")
  
  if (ploidy != 2 && dominance)
    stop("Dominance matrix only implemented for ploidy = 2.")
  
  if (slater && ploidy != 4) {
    warning(paste("Slater method is only valid for ploidy = 4.",
                  "Proceeding with default method."))
    slater <- FALSE
  }
  
  #-------------------------------------------------------------------------
  # 2. Pedigree Preprocessing and Verification
  #-------------------------------------------------------------------------
  orig.order <- as.character(data[[1]])
  
  if (verify) {
    message("Verifying conflicting data...\n")
    if (verifyped(data)) stop("Please double-check your data and try again.")
  }
  
  message("Organizing data with fast method...\n")
  data.after.treat <- try(datatreat(data = data, unk = 0, ...), silent = TRUE)
  
  if (inherits(data.after.treat, "try-error") ||
      !is.list(data.after.treat) ||
      is.null(data.after.treat$ind_data) ||
      length(unique(data.after.treat$ind_data)) != nrow(data)) {
    stop(paste("It wasn't possible to organize your data chronologically.",
               "Check for conflicting pedigree entries, duplicate IDs, or circular references."))
  }
  
  #-------------------------------------------------------------------------
  # 3. Matrix Construction (via Rcpp)
  #-------------------------------------------------------------------------
  data <- data.after.treat
  s <- as.integer(data$sire)
  d <- as.integer(data$dire)
  n <- length(s)
  
  if (!is.numeric(s) || !is.numeric(d) || length(s) != length(d))
    stop("Sire and dam columns must be numeric vectors of equal length.")
  
  if (n > 1000)
    message("Processing large pedigree data. It may take a couple of minutes...\n")
  
  start_time <- Sys.time()
  
  if (ploidy == 2) {
    message("Constructing matrix A using ploidy = 2\n")
    A <- buildA_ploidy2_cpp(s, d, n)
    if (dominance) {
      message("Constructing dominance relationship matrix\n")
      A <- buildDominanceMatrix_cpp(A, s, d)
    }
    
  } else if (slater && ploidy == 4) {
    message(
      sprintf(
        "Constructing matrix A using Slater et al. (2014), ploidy = 4, w = %.2f\n", 
        w))
    A <- buildA_slater_cpp(s, d, w)
    
  } else {
    message(
      sprintf(
        "Constructing matrix A using Kerr et al. (2012), ploidy = %d, w = %.2f\n", 
        ploidy, w))
    v <- ploidy / 2
    A <- buildA_kerr_cpp(s, d, w, v)
  }
  
  #-------------------------------------------------------------------------
  # 4. Post-processing
  #-------------------------------------------------------------------------
  if (anyNA(A))
    message(paste("Warning: Matrix contains NA values. Use",
                  "'verifyped()' to check data.\n"))
  
  rownames(A) <- colnames(A) <- data$ind.data
  A <- A[orig.order, orig.order]
  
  if (ASV) {
    A <- get_ASV(A)
  }
  
  attr(A, "ploidy") <- ploidy
  attr(A, "method") <- if (dominance) {
    "dominance"
  } else if (slater) {
    "slater"
  } else {
    "default"
  }
  
  elapsed <- difftime(Sys.time(), start_time, units = "mins")
  message("Completed! Time =", round(elapsed, 2), "minutes\n")
  
  return(A)
}
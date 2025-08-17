########################################
# 									
# Package: AGHmatrix 							
# 									
# File: Hmatrix.R 							
# Contains: Hmatrix 							
# 									
# Written by Rodrigo Rampazo Amadeu 					
# 									
# First version: Feb-2014 						
# Last update: 09-Jul-2019 						
# License: GPL-3
# 									
#######################################

#' Construction of Combined Relationship Matrix H
#'
#' Given a matrix A and a matrix G returns a H matrix. H matrix is the relationship matrix using combined information from the pedigree and genomic relationship matrices. First, you need to compute the matrices separated and then use them as input to build the combined H matrix. 
#' Two methods are implemented: `Munoz` shrinks the G matrix towards the A matrix scaling the molecular relatadness by each relationship classes; 
#' `Martini` is a modified version from Legarra et al. (2009) where combines A and G matrix using scaling factors. When method is equal `Martini` and `tau=1` and `omega=1` you have the same H matrix as in Legarra et al. (2009).
#'
#' @param A A matrix from function Amatrix
#' @param G G matrix from function Gmatrix
#' @param markers matrix marker which generated the Gmatrix
#' @param c constant value of H computation, default: c=0
#' @param method "Martini" or "Munoz", default="Martini"
#' @param missingValue missing value in data, default=-9.
#' @param maf max of missing data accepted to each markerm default=0.05.
#' @param ploidy data ploidy (an even number between 2 and 20), default=2.
#' @param tau to be used for Martini's method, default=1. 
#' @param omega to be used of Martini's method, default=1.
#' @param roundVar only used for Munoz's method, how many digits to consider the relationship be of same class, default=2.
#' @param ASV if TRUE, transform matrix into average semivariance (ASV) equivalent (K = K / (trace(K) / (nrow(K)-1))). Details formula 2 of Fieldmann et al. (2022). Default = FALSE.
#' 
#' @return H Matrix with the relationship between the individuals based on pedigree and corrected by molecular information
#'
#' @examples 
#' \dontrun{
#' data(ped.sol)
#' data(snp.sol)
#' #Computing the numerator relationship matrix 10% of double-reduction
#' Amat <- Amatrix(ped.sol, ploidy=4, w = 0.1)
#' #Computing the additive relationship matrix based on VanRaden (modified)
#' Gmat <- Gmatrix(snp.sol, ploidy=4, 
#'                 maf=0.05, method="VanRaden")
#' Gmat <- round(Gmat,3) #to be easy to invert
#' 
#' #Computing H matrix (Martini)
#' Hmat_Martini <- Hmatrix(A=Amat, G=Gmat, method="Martini", 
#'                      ploidy=4, 
#'                      maf=0.05)
#'                      
#' #Computing H matrix (Munoz)
#' Hmat_Munoz <- Hmatrix(A=Amat, G=Gmat, markers = snp.sol, 
#'                       ploidy=4, method="Munoz",
#'                       roundVar=2,
#'                       maf=0.05)
#' }
#'
#' @author Rodrigo R Amadeu, \email{rramadeu@@gmail.com}
#' @author Thiago de Paula Oliveira \email{toliveira@abacusbio.com}
#' @references \emph{Feldmann MJ, et al. 2022. Average semivariance directly yields accurate estimates of the genomic variance in complex trait analyses. G3 (Bethesda), 12(6).}
#' @references \emph{Munoz, PR. 2014 Unraveling additive from nonadditive effects using genomic relationship matrices. Genetics 198, 1759-1768}
#' @references \emph{Martini, JW, et al. 2018 The effect of the H-1 scaling factors tau and omega on the structure of H in the single-step procedure. Genetics Selection Evolution 50(1), 16}
#' @references \emph{Legarra, A, et al. 2009 A relationship matrix including full pedigree and genomic information. Journal of Dairy Science 92, 4656–4663}
#' @export
Hmatrix <- function(A=NULL,
                    G=NULL,
                    markers=NULL,
                    c=0,
                    method="Martini",
                    tau=1,
                    omega=1,
                    missingValue=-9,
                    maf=0,
                    ploidy=2,
                    roundVar=3,
                    ASV=FALSE){
  #-----------------------------------------------------------------------------
  # pre-checks
  #-----------------------------------------------------------------------------
  if (is.null(A) || is.null(G))
    stop("Please provide both A and G.")
  
  valid_methods <- c("Martini", "Munoz")
  if (!is.character(method) || length(method) != 1L || 
      !(method %in% valid_methods)) {
    stop("Invalid `method`: ", sQuote(method),
         ". Allowed methods: ", paste(sQuote(valid_methods), collapse = ", "), 
         ".", call. = FALSE)
  }
  
  # Basic input checks
  if (!is.matrix(A) || !is.matrix(G))
    stop("A and G must be matrices.")
  
  if (is.null(rownames(A)) || is.null(colnames(A)) ||
      is.null(rownames(G)) || is.null(colnames(G)))
    stop("A and G must be square with dimnames (IDs).")
  
  if (!identical(rownames(A), colnames(A)) ||
      !identical(rownames(G), colnames(G)))
    stop("A and G must be square with matching row/colnames in each.")
  
  Aorig <- A; Gorig <- G
  An <- rownames(Aorig); Gn <- rownames(Gorig)
  
  Time <- proc.time()
  cat("Comparing the matrices... \n")
  #-----------------------------------------------------------------------------
  # Missing pattern
  #-----------------------------------------------------------------------------
  # in A but not in G
  idxGmiss <- is.na(match(An, Gn))
  # in G but not in A
  idxAmiss <- is.na(match(Gn, An))
  if (any(!idxAmiss)) {
    # order we’ll keep in G submatrix
    Gnhat <- Gn[!idxAmiss]             
  } else Gnhat <- Gn
  if (any(!idxGmiss)) {
    # order we’ll keep in A submatrix (Munoz method)
    Anhat <- An[!idxGmiss]
  } else Anhat <- An
  
  #-----------------------------------------------------------------------------
  # keep G on its available IDs
  #-----------------------------------------------------------------------------
  G <- Gorig[Gnhat, Gnhat, drop = FALSE]
  elapsed <- (proc.time() - Time)[["elapsed"]]
  cat("Completed! Time =", elapsed / 60, "minutes\n")
  
  #-----------------------------------------------------------------------------
  # main dispatch
  #-----------------------------------------------------------------------------
  cat("Computing the H matrix... \n")
  Time <- proc.time()
  
  #-----------------------------------------------------------------------------
  # Martini
  #-----------------------------------------------------------------------------
  if (method == "Martini") {
    # build the global ID order
    idH <- unique(c(rownames(G), rownames(Aorig)))
    idH <- rev(idH)
    
    # reindex A once in that order
    A <- Aorig[idH, idH, drop = FALSE]
    inG <- !is.na(match(idH, rownames(G)))
    # partition indices
    idxAonly <- which(!inG)  # in A only
    idxG     <- which(inG)   # in G (genotyped)
    
    # pull the blocks we need
    A12 <- A[idxAonly, idxG, drop=FALSE]
    A22 <- A[idxG,     idxG, drop=FALSE]
    # same order as A22
    G22 <- G[idH[idxG], idH[idxG], drop=FALSE]
    
    # C++ core (does all inversions and multiplications)
    blk <- H_martini_blocks(A12, A22, G22, tau = tau, omega = omega)
    H11 <- blk$H11
    H12 <- blk$H12
    H21 <- blk$H21
    H22corr <- blk$H22corr
    
    # start from A and add corrections in-place by blocks
    H <- A
    if (length(idxAonly)) {
      H[idxAonly, idxAonly] <- H[idxAonly, idxAonly] + H11
      H[idxAonly, idxG]     <- H[idxAonly, idxG]     + H12
      H[idxG,     idxAonly] <- H[idxG,     idxAonly] + H21
    }
    if (length(idxG)) {
      H[idxG, idxG] <- H[idxG, idxG] + H22corr
    }
    
    if (ASV) H <- get_ASV(H)
    
    # attributes
    nm <- attr(Gorig, "nmarkers", exact = TRUE)
    if (is.null(nm)) nm <- NA_integer_
    attr(H, "method")   <- method
    attr(H, "ploidy")   <- ploidy
    attr(H, "nmarkers") <- nm
    
    Time <- as.matrix(proc.time()-Time)
    cat("\n","Completed! Time =", Time[3]/60," minutes \n")
    return(H)
  }
  
  #-----------------------------------------------------------------------------
  # Munoz
  #-----------------------------------------------------------------------------
  if (method == "Munoz") {
    # A in A-order on Anhat; G in G-order on Gnhat.
    A <- Aorig[Anhat, Anhat, drop = FALSE]
    G <- Gorig[Gnhat, Gnhat, drop = FALSE]
    
    # MarkersMatrix computed on raw markers, then subset in G’s order
    if (is.null(markers))
      stop("Aborting: For Munoz method you need to specify 'markers'.")
    Mm <- as.matrix(data.matrix(markers))
    if (!is.na(missingValue)) Mm[Mm == missingValue] <- NA
    markersmatrix <- Gmatrix_MarkersMask(Mm)
    # give it names (needed for character subsetting)
    if (!is.null(rownames(Mm))) {
      dimnames(markersmatrix) <- list(rownames(Mm), rownames(Mm))
    }
    if (!is.null(rownames(markersmatrix))) {
      markersmatrix <- markersmatrix[Gnhat, Gnhat, drop = FALSE]
    } else {
      if (nrow(markersmatrix) != length(Gnhat))
        stop("MarkersMatrix has no rownames and size doesn’t match G.")
    }
    
    # C++ calculation
    agg <- munoz_var_mean_by_Aclass(G, A, round_digits = roundVar)
    varA     <- agg$varA
    classes  <- as.numeric(agg$classes)
    v_by_cls <- as.numeric(agg$var_by_class)
    v_by_cls[v_by_cls == 0]         <- NA_real_
    v_by_cls[is.infinite(v_by_cls)] <- NA_real_
    if (anyNA(v_by_cls)) {
      filled <- approx(x = seq_along(v_by_cls), y = v_by_cls,
                       xout = seq_along(v_by_cls), rule = 2)$y
      # rebuild varA
      class_id <- match(round(A, roundVar), classes)
      varA <- matrix(filled[class_id], nrow(A), ncol(A),
                     dimnames = dimnames(A))
    }
    
    # beta and H
    if (any(markersmatrix <= 0, na.rm = TRUE))
      stop("MarkersMatrix contains zeros or negatives; cannot compute 1/count.")
    invM <- 1 / markersmatrix
    beta <- 1 - (c + invM / varA)
    Hloc <- beta * (G - A) + A
    
    # write back into A’s frame at Anhat (legacy)
    H <- Aorig
    H[Anhat, Anhat] <- Hloc
    
    if (ASV) H <- get_ASV(H)
    
    # attributes
    nm <- tryCatch(ncol(markers), error = function(e) NA_integer_)
    attr(H, "method")   <- method
    attr(H, "ploidy")   <- ploidy
    attr(H, "nmarkers") <- nm
    
    Time <- as.matrix(proc.time() - Time)
    cat("\n","Completed! Time =", Time[3]/60," minutes \n")
    cat("\n","Returning H = A matrix corrected by G... \n")
    return(H)
  }
}
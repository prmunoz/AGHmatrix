# Parity tests between legacy R maths and the new C++ 

# Legacy Helpers
Mcheck_l = function(SNPmatrix = NULL,
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
    mono = (apply(SNPmatrix, 2, var, na.rm=TRUE)==0)
    mono[is.na(mono)] = TRUE
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

get_ASV_l <- function(x){ x / ( sum(diag(x)) / (nrow(x) - 1)) }

slater_par_l <- function(X, ploidy){
  prime.index <- c(3,5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61,67,71,73,79)
  NumberMarkers <- ncol(X)
  nindTotal <- nrow(X)
  X <- X + 1
  temp <- seq(1, NumberMarkers, 10000)
  temp <- cbind(temp, temp + 9999)
  temp[length(temp)] <- NumberMarkers
  prime.index <- prime.index[1:(ploidy+1)]
  for(i in 1:nrow(temp)){
    X.temp <- X[, c(temp[i,1]:temp[i,2]), drop=FALSE]
    NumberMarkers <- ncol(X.temp)
    X.temp <- X.temp %*% t(kronecker(diag(NumberMarkers), prime.index))
    X.temp[which(as.vector(X.temp) %in% c(prime.index * c(1:(ploidy+1))))] <- 1
    X.temp[X.temp != 1] <- 0
    if(i==1){ X_out <- X.temp } else { X_out <- cbind(X_out, X.temp) }
  }
  gc()
  X_out
}

check_Gmatrix_data <- function(SNPmatrix, ploidy, method, ratio=FALSE, integer=TRUE){
  if (is.null(SNPmatrix)) stop("Please define the variable SNPdata")
  if (all(method != c("Yang", "VanRaden", "Slater", "Su", "Vitezica", "MarkersMatrix","Endelman"))) {
    stop("Unsupported method")
  }
  if(!is.matrix(SNPmatrix)) stop("SNPmatrix class must be matrix. Please verify it.")
  if(!ratio){
    if( ploidy > 20 | (ploidy %% 2) != 0) stop("Only even ploidy from 2 to 20")
    t <- max(SNPmatrix, na.rm = TRUE); if(t > ploidy) stop("Values above ploidy")
    t <- min(SNPmatrix, na.rm = TRUE); if(t < 0) stop("Values under 0")
    if(integer) if(prod(SNPmatrix == round(SNPmatrix), na.rm = TRUE)==0) stop("Non-integer values found")
  } else {
    t <- max(SNPmatrix, na.rm = TRUE); if(t > 1) stop("Ratio values above 1")
    t <- min(SNPmatrix, na.rm = TRUE); if(t < 0) stop("Ratio values under 0")
  }
}

Gmatrix_legacy <- function (SNPmatrix = NULL, method = "VanRaden", 
                            missingValue = -9, maf = 0, thresh.missing = .50,
                            verify.posdef = FALSE, ploidy=2,
                            pseudo.diploid = FALSE, integer=TRUE,
                            ratio = FALSE, impute.method = "mean", rmv.mono=FALSE, thresh.htzy=0,
                            ratio.check = TRUE, weights = NULL, ploidy.correction = FALSE, ASV=FALSE){
  Time <- proc.time(); markers <- colnames(SNPmatrix)
  if(!is.null(weights)) if(length(weights)!=ncol(SNPmatrix)) stop("weight should match #markers")
  if(ratio) method <- "VanRaden"
  if (!is.na(missingValue)) { m <- match(SNPmatrix, missingValue, 0); SNPmatrix[m > 0] <- NA }
  check_Gmatrix_data(SNPmatrix=SNPmatrix, method=method, ploidy=ploidy, ratio=ratio, integer=integer)
  NumberMarkers <- ncol(SNPmatrix)
  if(ratio==FALSE){
    SNPmatrix <- Mcheck_l(SNPmatrix, ploidy=ploidy, thresh.maf=maf, rmv.mono=rmv.mono,
                        thresh.htzy=thresh.htzy, thresh.missing=thresh.missing, impute.method=impute.method)
  }
  if(ratio && ratio.check){
    SNPmatrix <- Mcheck_l(SNPmatrix, ploidy=ploidy, thresh.maf=maf, rmv.mono=rmv.mono,
                        thresh.htzy=thresh.htzy, thresh.missing=thresh.missing, impute.method=impute.method)
  }
  if(method=="Slater"){
    P <- colSums(SNPmatrix,na.rm=TRUE)/nrow(SNPmatrix)
    SNPmatrix[, which(P>ploidy/2)] <- ploidy - SNPmatrix[, which(P>(ploidy/2))]
    SNPmatrix <- slater_par_l(SNPmatrix, ploidy=ploidy)
    NumberMarkers <- ncol(SNPmatrix)
    Frequency <- colSums(SNPmatrix, na.rm=TRUE)/nrow(SNPmatrix)
    FreqP <- matrix(rep(Frequency, each=nrow(SNPmatrix)), ncol=ncol(SNPmatrix))
  }
  if(ploidy==2){
    alelleFreq <- function(x, y) {(2 * length(which(x == y)) + length(which(x == 1)))/(2 * length(which(!is.na(x))))}
    Frequency <- cbind(apply(SNPmatrix, 2, function(x) alelleFreq(x,0)), apply(SNPmatrix, 2, function(x) alelleFreq(x,2)))
    FreqP <- matrix(rep(Frequency[, 2], each = nrow(SNPmatrix)), ncol = ncol(SNPmatrix))
  }
  if(ploidy>2 && pseudo.diploid){
    P <- colSums(SNPmatrix,na.rm = TRUE)/nrow(SNPmatrix)
    SNPmatrix[,which(P>ploidy/2)] <- ploidy-SNPmatrix[,which(P>(ploidy/2))]
    Frequency <- colSums(SNPmatrix,na.rm=TRUE)/(ploidy*nrow(SNPmatrix))
    Frequency <- cbind(1-Frequency,Frequency)
    FreqP <- matrix(rep(Frequency[, 2], each = nrow(SNPmatrix)), ncol = ncol(SNPmatrix))
    SNPmatrix[SNPmatrix %in% c(1:(ploidy-1))] <- 1
    SNPmatrix[SNPmatrix==ploidy] <- 2
  }
  if (method == "MarkersMatrix"){
    Gmatrix <- !is.na(SNPmatrix)
    Gmatrix <- tcrossprod(Gmatrix, Gmatrix)
    return(Gmatrix)
  }
  if (method == "VanRaden"){
    if(is.null(weights)){
      if(ploidy==2 & ratio==FALSE){
        TwoPQ <- 2 * t(Frequency[, 1]) %*% Frequency[, 2]
        SNPmatrix <- SNPmatrix - 2 * FreqP
        SNPmatrix[is.na(SNPmatrix)] <- 0
        Gmatrix <- (tcrossprod(SNPmatrix, SNPmatrix))/as.numeric(TwoPQ)
      } else {
        if(ploidy.correction){
          if(ratio==FALSE){
            Frequency <- apply(X=SNPmatrix, FUN=mean, MARGIN=2, na.rm=TRUE)/ploidy
            K <- sum(ploidy * Frequency * (1-Frequency))
          } else {
            Frequency <- apply(X=SNPmatrix, FUN=mean, MARGIN=2, na.rm=TRUE)
            K <- sum(1/ploidy * Frequency * (1-Frequency))
          }
        }
        SNPmatrix <- scale(SNPmatrix, center=TRUE, scale=FALSE)
        if(!ploidy.correction){ K <- sum(apply(X=SNPmatrix, FUN=var, MARGIN=2, na.rm=TRUE)) }
        SNPmatrix[which(is.na(SNPmatrix))] <- 0
        Gmatrix <- tcrossprod(SNPmatrix)/K
      }
    } else {
      weights <- weights[match(colnames(SNPmatrix), markers)]
      if(ploidy==2 & ratio==FALSE){
        TwoPQ <- 2 * t(Frequency[, 1]) %*% Frequency[, 2]
        SNPmatrix <- SNPmatrix - 2 * FreqP
        SNPmatrix[is.na(SNPmatrix)] <- 0
        Gmatrix <- tcrossprod(tcrossprod(SNPmatrix, diag(weights)), SNPmatrix)/as.numeric(TwoPQ)
      } else {
        if(ploidy.correction){
          if(ratio==FALSE){
            Frequency <- apply(X=SNPmatrix, FUN=mean, MARGIN=2, na.rm=TRUE)/ploidy
            K <- sum(ploidy * Frequency * (1-Frequency))
          } else {
            Frequency <- apply(X=SNPmatrix, FUN=mean, MARGIN=2, na.rm=TRUE)
            K <- sum(Frequency * (1-Frequency))
          }
        }
        SNPmatrix <- scale(SNPmatrix, center=TRUE, scale=FALSE)
        if(!ploidy.correction){ K <- sum(apply(X=SNPmatrix, FUN=var, MARGIN=2, na.rm=TRUE)) }
        SNPmatrix[which(is.na(SNPmatrix))] <- 0
        Gmatrix <- tcrossprod(tcrossprod(SNPmatrix, diag(weights)), SNPmatrix)/K
      }
    }
  }
  if (method == "Yang"){
    FreqPQ <- matrix(rep(2 * Frequency[, 1] * Frequency[, 2], each = nrow(SNPmatrix)), ncol = ncol(SNPmatrix))
    G.all <- (SNPmatrix^2 - (1 + 2 * FreqP) * SNPmatrix + 2 * (FreqP^2))/FreqPQ
    G.ii <- as.matrix(colSums(t(G.all), na.rm = TRUE))
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
    Gmatrix <- tcrossprod(SNPmatrix,SNPmatrix)/sum(TwoPQ[1,]*(1-TwoPQ[1,]))
  }
  if (method == "Vitezica"){
    TwoPQ <- 2*(FreqP[1,])*(1-FreqP[1,])
    SNPmatrix[is.na(SNPmatrix)] <- 0
    SNPmatrix <- (SNPmatrix==0)*-2*(FreqP^2) + (SNPmatrix==1)*2*(FreqP)*(1-FreqP) + (SNPmatrix==2)*-2*((1-FreqP)^2)
    Gmatrix <- tcrossprod(SNPmatrix,SNPmatrix)/sum(TwoPQ^2)
  }
  if (method == "Slater"){
    drop.alleles <- which(Frequency==0)
    if(length(drop.alleles)>0){
      Frequency <- Frequency[-drop.alleles]
      SNPmatrix <- SNPmatrix[,-drop.alleles, drop=FALSE]
      FreqP <- FreqP[,-drop.alleles, drop=FALSE]
    }
    FreqPQ <- matrix(rep(Frequency * (1-Frequency), each = nrow(SNPmatrix)), ncol = ncol(SNPmatrix))
    SNPmatrix[which(is.na(SNPmatrix))] <- 0
    G.ii <- (SNPmatrix^2 - (2 * FreqP) * SNPmatrix + FreqP^2)/FreqPQ
    G.ii <- as.matrix(colSums(t(G.ii), na.rm = TRUE))
    G.ii <- 1 + (G.ii)/NumberMarkers
    SNPmatrix <- (SNPmatrix - (FreqP))/sqrt(FreqPQ)
    SNPmatrix[is.na(SNPmatrix)] <- 0
    Gmatrix <- (tcrossprod(SNPmatrix, SNPmatrix))/NumberMarkers
    diag(Gmatrix) <- G.ii
  }
  if( method == "Endelman" ){
    if( ploidy != 4 ) stop("'Endelman' requires ploidy=4")
    Frequency <- colSums(SNPmatrix)/(nrow(SNPmatrix)*ploidy)
    Frequency <- cbind(Frequency,1-Frequency)
    SixPQ <- 6 * t((Frequency[, 1]^2)) %*% (Frequency[, 2]^2)
    SNPmatrix <- 6 * t((Frequency[, 1]^2)%*%t(rep(1,nrow(SNPmatrix)))) - 
      3*t((Frequency[, 1])%*%t(rep(1,nrow(SNPmatrix))))*SNPmatrix + 0.5 * SNPmatrix*(SNPmatrix-1)
    Gmatrix <- (tcrossprod(SNPmatrix, SNPmatrix))/as.numeric(SixPQ)
  }
  if(ASV){ Gmatrix <- get_ASV_l(Gmatrix) }
  Gmatrix
}

get_target_fun <- function(){ getOption("AGHmatrix_target_fun", Gmatrix) }

# ---------- Deterministic small data generators ----------
make_diploid <- function(){
  X <- cbind(
    c(0,1,2,0,1,2),     # M1
    c(0,0,1,1,2,2),     # M2
    c(2,1,0,2,1,0),     # M3
    c(1,1,1,1,1,1),     # M4 (monomorphic)
    c(0,-9,2,1,0,2),    # M5 (with missingValue=-9)
    c(2,2,1,1,0,0)      # M6
  )
  colnames(X) <- paste0("M", 1:ncol(X)); rownames(X) <- paste0("I", 1:nrow(X))
  storage.mode(X) <- "double"
  X
}

make_tetraploid <- function(){
  X <- cbind(
    c(0,1,2,3,4,2),   # T1
    c(4,3,2,1,0,2),   # T2
    c(0,0,2,2,4,2),   # T3 (even dosages)
    c(1,2,3,1,2,3),   # T4
    c(0,4,0,4,0,4),   # T5 (polarised)
    c(2,2,2,2,2,2)    # T6 (monomorphic)
  )
  colnames(X) <- paste0("T", 1:ncol(X)); rownames(X) <- paste0("J", 1:nrow(X))
  storage.mode(X) <- "double"
  X
}

make_ratio <- function(){
  X <- cbind(
    c(0.0, 0.2, 0.4, 0.6, 0.8),
    c(0.1, 0.1, 0.1, 0.9, 0.9),
    c(0.5, 0.5, 0.5, 0.5, 0.5),
    c(0.0, 0.3, 0.6, 0.3, 0.0)
  )
  colnames(X) <- paste0("R", 1:ncol(X)); rownames(X) <- paste0("K", 1:nrow(X))
  storage.mode(X) <- "double"
  X
}

strip_attrs <- function(x){
  at <- attributes(x)
  if (is.null(at)) return(x)
  keep <- c("dim","dimnames","class")
  attributes(x) <- at[names(at) %in% keep]
  x
}
compare_G <- function(X, args, method, tol=1e-8){
  legacy <- do.call(Gmatrix_legacy, c(list(SNPmatrix=X, method=method), args))
  target_fun <- get_target_fun()
  target <- do.call(target_fun, c(list(SNPmatrix=X, method=method), args))
  # drop non-structural attributes so we compare pure numerics
  legacy <- strip_attrs(legacy)
  target <- strip_attrs(target)
  expect_equal(target, legacy, tolerance = tol, scale = 1)
  expect_true(isSymmetric(target))
  expect_identical(dim(target), c(nrow(X), nrow(X)))
  invisible(list(legacy=legacy, target=target))
}

# Tests per method

test_that("VanRaden (diploid) parity",{
  X <- make_diploid()
  args <- list(missingValue=-9, ploidy=2, ratio=FALSE, integer=TRUE,
               rmv.mono=FALSE, thresh.htzy=0, thresh.missing=1, verify.posdef=FALSE)
  compare_G(X, args, method="VanRaden")
})

test_that("Yang (diploid) parity",{
  X <- make_diploid()
  args <- list(missingValue=-9, ploidy=2, ratio=FALSE, integer=TRUE,
               rmv.mono=FALSE, thresh.htzy=0, thresh.missing = 1,
               impute.method="mean")
  compare_G(X, args, method="Yang")
})

test_that("Su (diploid dominance) parity",{
  X <- make_diploid()
  args <- list(missingValue=-9, ploidy=2, ratio=FALSE, integer=TRUE,
               rmv.mono=FALSE, thresh.htzy=0, thresh.missing=1)
  compare_G(X, args, method="Su")
})

test_that("Vitezica (diploid dominance) parity",{
  X <- make_diploid()
  args <- list(missingValue=-9, ploidy=2, ratio=FALSE, integer=TRUE,
               rmv.mono=FALSE, thresh.htzy=0, thresh.missing=1)
  compare_G(X, args, method="Vitezica")
})

test_that("MarkersMatrix parity",{
  X <- make_diploid()
  args <- list(missingValue=-9, ploidy=2)
  compare_G(X, args, method="MarkersMatrix", tol=0) # exact integer equality expected
})

test_that("VanRaden (tetraploid) parity, ploidy.correction off/on",{
  X <- make_tetraploid()
  base_args <- list(missingValue=-9, ploidy=4, ratio=FALSE, integer=TRUE,
                    rmv.mono=FALSE, thresh.htzy=0, thresh.missing=1)
  compare_G(X, base_args, method="VanRaden")
  compare_G(X, c(base_args, list(ploidy.correction=TRUE)), method="VanRaden")
})

test_that("Slater (autotetraploid full model) parity",{
  X <- make_tetraploid()
  args <- list(missingValue=-9, ploidy=4, ratio=FALSE, integer=TRUE,
               rmv.mono=FALSE, thresh.htzy=0, thresh.missing=1)
  compare_G(X, args, method="Slater")
})

test_that("Endelman (autotetraploid digentic dominance) parity",{
  X <- make_tetraploid()
  args <- list(missingValue=-9, ploidy=4, ratio=FALSE, integer=TRUE,
               rmv.mono=FALSE, thresh.htzy=0, thresh.missing=1)
  compare_G(X, args, method="Endelman")
})

test_that("Pseudo-diploid parametrisation parity (ploidy=4)",{
  X <- make_tetraploid()
  args <- list(missingValue=-9, ploidy=4, pseudo.diploid=TRUE, ratio=FALSE, integer=TRUE,
               rmv.mono=FALSE, thresh.htzy=0, thresh.missing=1)
  compare_G(X, args, method="VanRaden")
})

test_that("Weights alignment remains correct after marker filtering (no filtering here)",{
  X <- make_diploid()
  w <- c(M1=0.5, M2=1.2, M3=0.7, M4=1.0, M5=0.8, M6=1.5)
  # permute column order to ensure matching is by name, not position
  X2 <- X[, c("M3","M1","M6","M2","M5","M4")]
  args <- list(missingValue=-9, ploidy=2, ratio=FALSE, integer=TRUE,
               rmv.mono=FALSE, thresh.htzy=0, thresh.missing=1, weights=w)
  compare_G(X2, args, method="VanRaden")
})

test_that("ASV rescaling parity wraps any method",{
  X <- make_diploid()
  args <- list(missingValue=-9, ploidy=2, ratio=FALSE, integer=TRUE,
               rmv.mono=FALSE, thresh.htzy=0, thresh.missing=1, ASV=TRUE)
  compare_G(X, args, method="Yang")
  compare_G(X, args, method="VanRaden")
})

test_that("Ratio pathway forces VanRaden and matches explicit VanRaden+ratio",{
  X <- make_ratio()
  args1 <- list(missingValue=-9, ploidy=2, ratio=TRUE, ratio.check=FALSE, integer=FALSE)
  # method argument ignored when ratio=TRUE; we still pass Yang to check override
  out1 <- compare_G(X, args1, method="Yang")
  args2 <- list(missingValue=-9, ploidy=2, ratio=TRUE, ratio.check=FALSE, integer=FALSE)
  out2 <- compare_G(X, args2, method="VanRaden")
  expect_equal(out1$target, out2$target, tolerance=0)
  expect_equal(out1$legacy, out2$legacy, tolerance=0)
})

prop_methods <- c("VanRaden","Yang","Su","Vitezica")

test_that("All diploid methods return symmetric matrices with correct size",{
  X <- make_diploid(); args <- list(missingValue=-9, ploidy=2, ratio=FALSE, integer=TRUE,
                                    rmv.mono=FALSE, thresh.htzy=0, thresh.missing=1)
  for(m in prop_methods){
    target <- do.call(get_target_fun(), c(list(SNPmatrix=X, method=m), args))
    expect_true(isSymmetric(target)); expect_identical(dim(target), c(nrow(X), nrow(X)))
  }
})

test_that("new Gmatrix matches legacy on snp.pine", {
  # Load data
  data(snp.pine, package = "AGHmatrix")
  
  # Do a fake imputation for the second comparison
  snp.pine.fake <- ifelse(snp.pine == -9, 1, snp.pine)
  
  # New implementation
  G_VanRadenPine  <- Gmatrix(SNPmatrix = snp.pine,      missingValue = -9,
                             maf = 0, method = "VanRaden")
  G_VanRadenPine2 <- Gmatrix(SNPmatrix = snp.pine.fake, missingValue = -9,
                             maf = 0, method = "VanRaden")
  
  # Legacy implementation
  G_VanRadenPine_curr  <- Gmatrix_legacy(SNPmatrix = snp.pine,      missingValue = -9,
                                         maf = 0, method = "VanRaden")
  G_VanRadenPine2_curr <- Gmatrix_legacy(SNPmatrix = snp.pine.fake, missingValue = -9,
                                         maf = 0, method = "VanRaden")
  
  expect_false(any(round(G_VanRadenPine_curr,  4) != round(G_VanRadenPine,  4)),
               info = {
                 d <- round(G_VanRadenPine_curr,4) - round(G_VanRadenPine,4)
                 paste0("Rounded(4) mismatch in pine with missings. ",
                        "max|delta|=", max(abs(d)), 
                        ", nnz(delta)=", sum(d != 0))
               })
  
  expect_false(any(round(G_VanRadenPine2_curr, 4) != round(G_VanRadenPine2, 4)),
               info = {
                 d <- round(G_VanRadenPine2_curr,4) - round(G_VanRadenPine2,4)
                 paste0("Rounded(4) mismatch in fake-imputed pine. ",
                        "max|delta|=", max(abs(d)), 
                        ", nnz(delta)=", sum(d != 0))
               })
  
  # symmetry checks ---
  expect_true(isSymmetric(G_VanRadenPine))
  expect_true(isSymmetric(G_VanRadenPine2))
  expect_true(isSymmetric(G_VanRadenPine_curr))
  expect_true(isSymmetric(G_VanRadenPine2_curr))
})

test_that("new Gmatrix matches legacy (Su) on snp.pine", {
  # Load data
  data(snp.pine, package = "AGHmatrix")
  
  # Do a fake imputation for the second comparison
  snp.pine.fake <- ifelse(snp.pine == -9, 1, snp.pine)
  
  # New implementation
  G_SuPine  <- Gmatrix(SNPmatrix = snp.pine,      missingValue = -9,
                       maf = 0.05, method = "Su")
  G_SuPine2 <- Gmatrix(SNPmatrix = snp.pine.fake, missingValue = -9,
                       maf = 0.05, method = "Su")
  
  # Legacy implementation
  G_SuPine_curr  <- Gmatrix_legacy(SNPmatrix = snp.pine,      missingValue = -9,
                                   maf = 0.05, method = "Su")
  G_SuPine2_curr <- Gmatrix_legacy(SNPmatrix = snp.pine.fake, missingValue = -9,
                                   maf = 0.05, method = "Su")
  
  # Equality up to rounding(4)
  expect_false(any(round(G_SuPine_curr,  4) != round(G_SuPine,  4)),
               info = {
                 d <- round(G_SuPine_curr,4) - round(G_SuPine,4)
                 paste0("Rounded(4) mismatch in pine (Su) with missings. ",
                        "max|delta|=", max(abs(d)), 
                        ", nnz(delta)=", sum(d != 0))
               })
  expect_false(any(round(G_SuPine2_curr, 4) != round(G_SuPine2, 4)),
               info = {
                 d <- round(G_SuPine2_curr,4) - round(G_SuPine2,4)
                 paste0("Rounded(4) mismatch in fake-imputed pine (Su). ",
                        "max|delta|=", max(abs(d)), 
                        ", nnz(delta)=", sum(d != 0))
               })
  
  # Symmetry checks
  expect_true(isSymmetric(G_SuPine))
  expect_true(isSymmetric(G_SuPine2))
  expect_true(isSymmetric(G_SuPine_curr))
  expect_true(isSymmetric(G_SuPine2_curr))
})


test_that("new Gmatrix matches legacy (Yang) on snp.pine", {
  # Load data
  data(snp.pine, package = "AGHmatrix")
  
  # Do a fake imputation for the second comparison
  snp.pine.fake <- ifelse(snp.pine == -9, 1, snp.pine)
  
  # New implementation
  G_YangPine  <- Gmatrix(SNPmatrix = snp.pine,      missingValue = -9,
                         maf = 0.05, method = "Yang")
  G_YangPine2 <- Gmatrix(SNPmatrix = snp.pine.fake, missingValue = -9,
                         maf = 0.05, method = "Yang")
  
  # Legacy implementation
  G_YangPine_curr  <- Gmatrix_legacy(SNPmatrix = snp.pine,      missingValue = -9,
                                     maf = 0.05, method = "Yang")
  G_YangPine2_curr <- Gmatrix_legacy(SNPmatrix = snp.pine.fake, missingValue = -9,
                                     maf = 0.05, method = "Yang")
  
  # Equality up to rounding(4)
  expect_false(any(round(G_YangPine_curr,  4) != round(G_YangPine,  4)),
               info = {
                 d <- round(G_YangPine_curr,4) - round(G_YangPine,4)
                 paste0("Rounded(4) mismatch in pine (Yang) with missings. ",
                        "max|delta|=", max(abs(d)), 
                        ", nnz(delta)=", sum(d != 0))
               })
  expect_false(any(round(G_YangPine2_curr, 4) != round(G_YangPine2, 4)),
               info = {
                 d <- round(G_YangPine2_curr,4) - round(G_YangPine2,4)
                 paste0("Rounded(4) mismatch in fake-imputed pine (Yang). ",
                        "max|delta|=", max(abs(d)), 
                        ", nnz(delta)=", sum(d != 0))
               })
  
  # Symmetry checks
  expect_true(isSymmetric(G_YangPine))
  expect_true(isSymmetric(G_YangPine2))
  expect_true(isSymmetric(G_YangPine_curr))
  expect_true(isSymmetric(G_YangPine2_curr))
})

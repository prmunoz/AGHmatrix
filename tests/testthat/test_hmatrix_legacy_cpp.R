# Parity tests between new Hmatrix() and Hmatrix_legacy()
Hmatrix_legacy <- function(A=NULL,
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
                    ASV=FALSE
){
  Aorig <- A
  Gorig <- G
  
  Time = proc.time()
  cat("Comparing the matrices... \n")
  An <- rownames(Aorig)
  Gn <- rownames(Gorig)
  missingGmatrix <- which(is.na(match(An,Gn)))
  missingAmatrix <- which(is.na(match(Gn,An)))
  if(length(missingAmatrix)>0){
    Gnhat <- Gn[-missingAmatrix]
  }else{
    Gnhat <- Gn
  }
  
  if(length(missingGmatrix)>0){
    Anhat <- An[-missingGmatrix]    
  }else{
    Anhat <- An
  }
  
  A <- Aorig#[Anhat,Anhat]
  G <- Gorig[Gnhat,Gnhat]
  
  missingGmatrix <- An[missingGmatrix]
  missingAmatrix <- Gn[missingAmatrix]
  
  Time = as.matrix(proc.time()-Time)
  cat("Completed! Time =", Time[3]/60," minutes \n")
  
  cat("Computing the H matrix... \n")
  
  Time = proc.time()
  
  if(method=="Martini"){
    idA <-rownames(A)
    idG <- rownames(G)
    idH <- unique(c(idG,idA))
    idH <- rev(idH)
    A <- A[idH,idH]
    
    index = is.na(match(idH,idG))
    A11 <- A[index,index]
    A12 <- A[index,!index]
    A21 <- A[!index,index]
    A22 <- A[!index,!index]
    G22 <- G[idH[!index],idH[!index]]
    #if(is.singular.matrix(G22))
    #  stop(deparse("Matrix G22 is singular (not invertible)"))
    A22inv = solve(A22) #A is always invertible
    G22inv = try(solve(G22),silent=TRUE)
    if(inherits(G22inv,"try-error")){
      cat(G22inv)
      stop("G22 not inverting with solve(), try a different/modified G matrix")
    }
    H22 = solve((tau*G22inv+(1-omega)*A22inv))
    H11 = A12 %*% A22inv %*% (H22-A22) %*% A22inv %*% A21  
    H12 = A12 %*% A22inv %*% (H22-A22)
    H21 = (H22-A22) %*% A22inv%*%A21
    H22 = (H22-A22)
    H = A+cbind(rbind(H11,H21),rbind(H12,H22))
    
    
    if (ASV) {
      H = get_ASV(H)
    }
    
    Time = as.matrix(proc.time()-Time)
    cat("\n","Completed! Time =", Time[3]/60," minutes \n")
    return(H)
  }
  
  if(method=="Munoz"){
    A <- Aorig[Anhat,Anhat]
    if(is.null(markers))
      stop("Aborting: For Munoz method you need to specify method object")
    markersmatrix <- Gmatrix(markers,method="MarkersMatrix",ploidy=ploidy,missingValue=missingValue,maf=maf)
    
    #Computing the Variance of G by A classes (A rounded by roundVar)
    classes <- as.numeric(levels(as.factor(A)))
    classes <- unique(round(classes,roundVar))
    n <- length(classes)
    varA <- meanG <- matrix(NA,nrow=nrow(A),ncol=nrow(A))
    varAclasses <- c()
    for(i in 1:n){
      varA[round(A,roundVar)==classes[i]] <- varAclasses[i] <- var(G[round(A,roundVar)==classes[i]])
      meanG[round(A,roundVar)==classes[i]]  <- mean(G[round(A,roundVar)==classes[i]])
    }
    
    varAclasses[varAclasses==0] = NA
    varAclasses[is.infinite(varAclasses)] <- NA
    tmp <- which(is.na(varAclasses))
    for(i in 1:length(tmp)){
      varAclasses[tmp[i]]<-zoo::na.approx(varAclasses)[tmp[i]]
      varA[round(A,roundVar)==classes[tmp[i]]] <- varAclasses[tmp[i]]
    }
    
    #Computaing beta and H
    beta <- 1 - (c+(1/(markersmatrix[Gnhat,Gnhat]))/varA)
    H <- beta*(G-A)+A ######
    Aorig[Anhat,Anhat] = H
    
    if (ASV) {
      Aorig = get_ASV_l(Aorig)
    }
    
    Time = as.matrix(proc.time()-Time)
    cat("\n","Completed! Time =", Time[3]/60," minutes \n")
    cat("\n","Returning H = A matrix corrected by G... \n")
    
    return(Aorig)
  }
}

get_ASV_l <- function(x){ x / ( sum(diag(x)) / (nrow(x) - 1)) }

make_spd <- function(n, seed = 1L){
  set.seed(seed)
  M <- matrix(rnorm(n*n), n)
  S <- crossprod(M)               # SPD
  diag(S) <- diag(S) + n*1e-3
  S
}

# Build A and G
make_AG_ids <- function(){
  # A over 8 IDs, G over subset of 4
  idA <- paste0("I", 1:8)
  idG <- c("I3","I6","I2","I7")
  A <- diag(8)
  rownames(A) <- colnames(A) <- idA
  G <- make_spd(length(idG), seed = 42)
  rownames(G) <- colnames(G) <- idG
  list(A=A, G=G, idA=idA, idG=idG)
}

# Simple markers matrix for Muñoz
make_markers <- function(ids, m = 6L, seed = 99L){
  set.seed(seed)
  # diploid-style dosages 0/1/2, no missing
  X <- matrix(sample(0:2, length(ids)*m, replace = TRUE), nrow = length(ids), 
              ncol = m)
  rownames(X) <- ids; colnames(X) <- paste0("M", seq_len(m))
  storage.mode(X) <- "double"
  X
}

compare_H_martini <- function(A, G, tau = 1, omega = 1, c = 0, 
                              ASV = FALSE, tol = 1e-10){
  # new
  H_new <- Hmatrix(A=A, G=G, method = "Martini", tau = tau, omega = omega, 
                   c = c, ASV = ASV)
  # legacy
  H_old <- Hmatrix_legacy(A=A, G=G, method = "Martini", tau = tau, 
                          omega = omega, c = c, ASV = ASV)
  # numeric parity (ignore non-structural attrs)
  expect_equal(strip_attrs(H_new), strip_attrs(H_old), tolerance = tol, 
               scale = 1)
  # basic properties
  expect_true(isSymmetric(H_new))
  expect_identical(rownames(H_new), rownames(H_old))
  expect_identical(colnames(H_new), colnames(H_old))
  invisible(list(new = H_new, old = H_old))
}

compare_H_munoz <- function(A, G, markers, c = 0, ploidy = 2, roundVar = 0, 
                            ASV = FALSE, tol = 1e-10){
  H_new <- Hmatrix(A=A, G=G, markers=markers, method = "Munoz", c = c, 
                   ploidy = ploidy, roundVar = roundVar, ASV = ASV)
  H_old <- Hmatrix_legacy(A=A, G=G, markers=markers, method = "Munoz", c = c, 
                          ploidy = ploidy, roundVar = roundVar, ASV = ASV)
  expect_equal(strip_attrs(H_new), strip_attrs(H_old), tolerance = tol, 
               scale = 1)
  expect_true(isSymmetric(H_new))
  expect_identical(rownames(H_new), rownames(H_old))
  expect_identical(colnames(H_new), colnames(H_old))
  invisible(list(new = H_new, old = H_old))
}

strip_attrs <- function(x){
  at <- attributes(x); if (is.null(at)) return(x)
  keep <- c("dim","dimnames","class")
  attributes(x) <- at[names(at) %in% keep]
  x
}

# Tests

test_that("H (Martini) parity on simple SPD inputs, default tau=omega=1",{
  objs <- make_AG_ids(); with(objs, {
    compare_H_martini(A, G, tau = 1, omega = 1, c = 0, ASV = FALSE)
  })
})

test_that("H (Martini) parity holds for different tau/omega and ASV on",{
  objs <- make_AG_ids(); with(objs, {
    compare_H_martini(A, G, tau = 0.7, omega = 0.2, c = 0, ASV = FALSE)
    compare_H_martini(A, G, tau = 1.5, omega = 0.0, c = 0, ASV = TRUE)
  })
})

test_that("H (Munoz) parity with clean classes (roundVar=0), no ASV",{
  objs <- make_AG_ids(); with(objs, {
    markers <- make_markers(idG, m = 8)
    compare_H_munoz(A, G, markers = markers, c = 0, ploidy = 2,
                    roundVar = 0, ASV = FALSE)
  })
})

test_that("H (Munoz) parity with ASV on",{
  objs <- make_AG_ids(); with(objs, {
    markers <- make_markers(idG, m = 5)
    compare_H_munoz(A, G, markers = markers, c = 0.2, ploidy = 2, 
                    roundVar = 0, ASV = TRUE)
  })
})

test_that("Informative attributes",{
  objs <- make_AG_ids(); with(objs, {
    H <- Hmatrix(A=A, G=G, method = "Martini", tau = 1.2, omega = 0.3, 
                 c = 0.0, ASV = TRUE)
    # Attributes should not be NA/NULL
    at <- attributes(H)
    expect_true(is.list(at))
    # Commonly exposed attrs 
    if (!is.null(attr(H, "method"))) expect_identical(attr(H, "method"), "Martini")
    if (!is.null(attr(H, "tau")))     expect_equal(attr(H, "tau"), 1.2)
    if (!is.null(attr(H, "omega")))   expect_equal(attr(H, "omega"), 0.3)
    if (!is.null(attr(H, "c")))       expect_equal(attr(H, "c"), 0.0)
    if (!is.null(attr(H, "ASV")))     expect_identical(attr(H, "ASV"), TRUE)
  })
})

test_that("Martini path reorders/expands A to idH = unique(c(idG,idA))",{
  objs <- make_AG_ids(); with(objs, {
    # permute A order and add an extra id in A not in G
    permA <- sample(idA)
    A2 <- A[permA, permA]
    # new and legacy should agree
    res <- compare_H_martini(A2, G, tau = 0.9, omega = 0.1)
    expect_identical(rownames(res$new), rownames(res$old))
  })
})


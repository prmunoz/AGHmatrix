# Parity tests between new Amatrix() and Amatrix_legacy()
Amatrix_legacy <- function(data = NULL,
                    ploidy=2,
                    w=0,
                    verify=TRUE,
                    dominance=FALSE,
                    slater=FALSE,
                    ASV=FALSE,
                    ...){
  if(ploidy%%2!=0)
    stop(deparse("Ploidy should be an even number"))
  
  if(ploidy!=2 & dominance)
    stop(deparse("Dominance relationship matrix is implemented only for ploidy=2"))
  
  if( is.null(data))
    stop(deparse("Please define the variable data"))
  
  unk=0
  cat("Verifying conflicting data... \n")
  flag<-verifyped_legacy(data)
  if(flag) 
    stop(deparse("Please double-check your data and try again"))
  
  cat("Organizing data... \n")
  orig.order <- as.character(data[,1])
  data.after.treat <- try(datatreat_legacy(data=data,unk=unk,...),silent=TRUE)
  
  # checking if order was fixed
  flag = FALSE
  flag = inherits(data.after.treat,"try-error")
  if(!flag)
    flag = (length(unique(data.after.treat$ind.data))!=nrow(data))
  
  if(flag){
    cat("To organize the data in a fast way wasn't possible... \n")
    cat("Trying to organize in a slow (naive) way... \n")
    data.sorted <- sortped_legacy(data)
    data.after.treat <- try(datatreat_legacy(data=data.sorted,unk=unk,...))
    
    # checking if order was fixed
    flag = FALSE
    flag = inherits(data.after.treat,"try-error")
    if(!flag)
      flag = (length(unique(data.after.treat$ind.data))!=nrow(data))
    
    if(flag){
      cat("It wasn't possible to organize your data chronologically. We recommend you to do it by hand or use the flag 'naive_sort=TRUE'. If the problem persists, please contact this package mainteiner \n")
      return()
    }
  }
  
  data <- data.after.treat
  
  s <- data$sire
  d <- data$dire
  
  if( is.null(s) || is.null(d) )
    stop(deparse("Please define the variable s (sire) and/or d (dire)"))
  
  if( length(s) != length(d) )
    stop(deparse("Please verify the variable s (sire) and/or d (dire), they don't have the same length"))
  
  if( !is.numeric(s) || !is.numeric(d) )
    stop(deparse("Pleasy verify your data, it has to be 2 numeric vectors"))
  
  if( length(data$sire) > 1000 )
    cat("Processing a large pedigree data... It may take a couple of minutes... \n")
  
  
  n <- length(s)
  A <- matrix(NA,ncol=n,nrow=n)
  
  Time = proc.time()
  
  #### For ploidy 2 ####
  if(ploidy == 2){
    w <- NA
    cat("Constructing matrix A using ploidy = 2 \n")
    A[1,1] <- 1
    for( i in 2:n){
      
      ## Both are unknown
      if( s[i] == 0 && d[i] == 0 ){
        A[i,i] <- 1
        for( j in 1:(i-1))
          A[j,i] <- A[i,j] <- 0
      }
      
      ## Sire is unknown
      if( s[i] == 0 && d[i] != 0 ){
        A[i,i] <- 1
        for( j in 1:(i-1))
          A[j,i] <- A[i,j] <- 0.5*(A[j,d[i]])
      }
      
      ## Dire is unknown
      if( d[i] == 0 && s[i] != 0 ){
        A[i,i] <- 1
        for( j in 1:(i-1))
          A[j,i] <- A[i,j] <- 0.5*(A[j,s[i]])
      }
      
      ## Both are known
      if( d[i] != 0 && s[i] != 0 ){
        A[i,i] <- 1+0.5*(A[d[i],s[i]])
        for( j in 1:(i-1))
          A[j,i] <- A[i,j] <- 0.5*(A[j,s[i]]+A[j,d[i]])
      }
    }
  }
  
  if(dominance){
    cat("Constructing dominance relationship matrix \n")
    D <- matrix(NA,ncol=n,nrow=n)
    for(i in 1:n){
      for(j in 1:n){
        u1 <- ifelse(length(A[s[i],s[j]])>0,A[s[i],s[j]],0)
        u2 <- ifelse(length(A[d[i],d[j]])>0,A[d[i],d[j]],0)
        u3 <- ifelse(length(A[s[i],d[j]])>0,A[s[i],d[j]],0)
        u4 <- ifelse(length(A[s[j],d[i]])>0,A[s[j],d[i]],0)
        D[i,j] <- D[j,i] <- 0.25*(u1*u2+u3*u4)
      }
    }
    diag(D)<-1
    A<-D
    D<-NULL
  }
  
  #### For ploidy 4 ####
  if(slater==TRUE){
    listA <- list()
    
    cat(paste("Constructing matrix A using ploidy = 4 and proportion of double reduction =",w,";as in Slater et al. (2014) \n"))
    start.time <- Sys.time()
    
    A[1,1] <- (1+w)/4
    for( i in 2:n){
      
      ## Both are unknown
      if( s[i] == 0 && d[i] == 0 ){
        A[i,i] <- (1+w)/4
        for( j in 1:(i-1))
          A[j,i] <- A[i,j] <- 0
      }
      
      ## Sire is unknown
      if( s[i] == 0 && d[i] != 0 ){
        A[i,i] <- (5 + 7*w + 4*A[d[i],d[i]]*(1-w) ) / 24
        for( j in 1:(i-1))
          A[j,i] <- A[i,j] <- 0.5*(A[j,d[i]])
      }
      
      ## Dire is unknown
      if( d[i] == 0 && s[i] != 0 ){
        A[i,i] <- (5 + 8*w + 4*A[s[i],s[i]]*(1-w) ) / 24 ##On Slater in 7w, deriving from hand based on Kerr is 8w
        for( j in 1:(i-1))
          A[j,i] <- A[i,j] <- 0.5*(A[j,s[i]])
      }
      
      ## Both are known
      if( d[i] != 0 && s[i] != 0 ){
        A[i,i] <- (1 + 2*w + (1-w)*(A[s[i],s[i]]) + (1-w)*(A[d[i],d[i]]) + 3*A[s[i],d[i]] ) / 6
        for( j in 1:(i-1))
          A[j,i] <- A[i,j] <- 0.5*(A[j,s[i]]+A[j,d[i]])
      }
    }
    
    A <- 4*A
  }
  
  if(slater==FALSE && ploidy>2){ ## It does not use double-reduction proportion, need to double-check formula on kerr 2012 for higher ploidies...
    listA <- list()
    cat(paste("Constructing matrix A using ploidy =",ploidy,"and proportion of double reduction =",w,";as in Kerr et al. (2012) \n"))
    start.time <- Sys.time()
    v = ploidy/2
    A[1,1] <- (1)/(2*v)
    for( i in 2:n){
      ## Both are unknown
      if( s[i] == 0 && d[i] == 0 ){
        A[i,i] <- (1)/(2*v)
        for( j in 1:(i-1))
          A[j,i] <- A[i,j] <- 0
      }
      
      ## Sire is unknown
      if( s[i] == 0 && d[i] != 0 ){
        A[i,i] <- (1 + (v-1)*w + ((v-1)*(1-w)*(v*A[d[i],d[i]] + 1/2 - 1))/(2*v-1))/(2*v)
        for( j in 1:(i-1))
          A[j,i] <- A[i,j] <- 0.5*(A[j,d[i]])
      }
      
      ## Dire is unknown
      if( d[i] == 0 && s[i] != 0 ){
        A[i,i] <- (1 + (v-1)*w + ((v-1)*(1-w)*(v*A[s[i],s[i]] + 1/2 - 1))/(2*v-1))/(2*v)
        for( j in 1:(i-1))
          A[j,i] <- A[i,j] <- 0.5*(A[j,s[i]])
      }
      
      ## Both are known
      if( d[i] != 0 && s[i] != 0 ){
        A[i,i] <- (1  + (v-1)*w + ((v-1)*(1-w)*(v*A[d[i],d[i]] + v*A[s[i],s[i]] - 1)/(2*v-1)))/(2*v) + A[d[i],s[i]]/2
        for( j in 1:(i-1))
          A[j,i] <- A[i,j] <- 0.5*(A[j,s[i]]+A[j,d[i]])
      }
    }
    
    A <- 2*v*A
  }
  
  NA.errors <- which(is.na(A))
  if( length(NA.errors) > 0 )
    cat("Please verify your original data with the function 'verifyped', there are some data missing/conflicting data \n")
  
  Time = as.matrix(proc.time()-Time)
  cat("Completed! Time =", Time[3]/60," minutes \n")
  #      "Visualization options: (matrix, w) \n ")
  rownames(A) <- colnames(A) <- data$ind.data
  
  A <- A[orig.order,orig.order]
  
  
  if (ASV) {
    A = get_ASV_legacy(A)
  }
  
  return(A)
}

# Helpers
verifyped_legacy <- function(data) FALSE
sortped_legacy <- function(data) data
datatreat_legacy <- function(data, unk = 0, ...) {
  stopifnot(is.data.frame(data), ncol(data) >= 3)
  colnames(data)[1:3] <- c("ind", "sire", "dam")
  ids <- as.character(data$ind)
  si  <- as.character(data$sire)
  da  <- as.character(data$dam)
  # Unknowns indicated by "0" or 0 -> keep as 0
  map <- setNames(seq_along(ids), ids)
  sire_idx <- ifelse(si == "0" | si == 0, 0L, as.integer(map[si]))
  dam_idx  <- ifelse(da == "0" | da == 0, 0L, as.integer(map[da]))
  out <- data.frame(
    ind.data = ids,
    sire     = sire_idx,
    dire     = dam_idx,
    stringsAsFactors = FALSE
  )
  out
}
get_ASV_legacy <- function(x) x / (sum(diag(x)) / (nrow(x) - 1))

# Pedigree
# Founders: A,B,C,D ; Progeny: E=A×B, F=C×D ; G=E×F ; H=E×F
make_ped <- function() {
  data.frame(
    id   = c("A","B","C","D","E","F","G","H"),
    sire = c("0","0","0","0","A","C","E","E"),
    dam  = c("0","0","0","0","B","D","F","F"),
    stringsAsFactors = FALSE
  )
}

make_ped2 <- function() {
  base <- make_ped()
  extra <- data.frame(
    id   = c("I","J"),
    sire = c("G","G"),
    dam  = c("H","E"),
    stringsAsFactors = FALSE
  )
  rbind(base, extra)
}

# ---- Utility: strip non-structural attrs and compare matrices ----
strip_attrs <- function(x) {
  at <- attributes(x)
  if (is.null(at)) return(x)
  keep <- c("dim", "dimnames", "class")
  attributes(x) <- at[names(at) %in% keep]
  x
}

compare_A <- function(ped, ..., tol = 1e-10) {
  legacy <- Amatrix_legacy(ped, ...)
  target <- Amatrix(ped, ...)
  legacy <- strip_attrs(legacy)
  target <- strip_attrs(target)
  expect_equal(target, legacy, tolerance = tol, scale = 1)
  expect_true(isSymmetric(target))
  expect_identical(dim(target), c(nrow(ped), nrow(ped)))
  invisible(list(legacy = legacy, target = target))
}

# Tests

test_that("Diploid additive (Henderson 1976) parity", {
  ped <- make_ped()
  compare_A(ped, ploidy = 2, dominance = FALSE, verify = FALSE)
})

test_that("Diploid dominance (Cockerham 1954) parity", {
  ped <- make_ped()
  compare_A(ped, ploidy = 2, dominance = TRUE, verify = FALSE)
})

test_that("ASV rescaling parity (diploid additive)", {
  ped <- make_ped()
  compare_A(ped, ploidy = 2, dominance = FALSE, ASV = TRUE, verify = FALSE)
})

test_that("Autotetraploid additive (Kerr 2012) parity, w = 0", {
  ped <- make_ped()
  compare_A(ped, ploidy = 4, w = 0, slater = FALSE, verify = FALSE)
})

test_that("Autotetraploid additive (Kerr 2012) parity, w > 0", {
  ped <- make_ped()
  compare_A(ped, ploidy = 4, w = 0.1, slater = FALSE, verify = FALSE)
})

test_that("Autotetraploid additive (Slater 2013) parity, w > 0", {
  ped <- make_ped()
  compare_A(ped, ploidy = 4, w = 0.1, slater = TRUE, verify = FALSE)
})

test_that("Higher ploidy (hexaploid, Kerr 2012) parity, w = 0", {
  ped <- make_ped2()
  compare_A(ped, ploidy = 6, w = 0, slater = FALSE, verify = FALSE)
})

test_that("Higher ploidy (hexaploid, Kerr 2012) parity, w > 0", {
  ped <- make_ped2()
  compare_A(ped, ploidy = 6, w = 0.05, slater = FALSE, verify = FALSE)
})

test_that("Row/column names and order preserved", {
  ped <- make_ped2()
  out1 <- Amatrix_legacy(ped, ploidy = 2, verify = FALSE)
  out2 <- Amatrix(ped, ploidy = 2, verify = FALSE)
  expect_identical(rownames(out2), rownames(out1))
  expect_identical(colnames(out2), colnames(out1))
})

test_that("Dominance requested with ploidy != 2", {
  ped <- make_ped()
  expect_error(Amatrix_legacy(ped, ploidy = 4, dominance = TRUE, verify = FALSE))
  expect_error(Amatrix(ped, ploidy = 4, dominance = TRUE, verify = FALSE))
})

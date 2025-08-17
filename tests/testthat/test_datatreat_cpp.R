datatreat_legacy <- function(data=NULL,
                      n.max=50,
                      unk=0,
                      save=FALSE
){
  indicator <- k <-  0
  if(is.null(data))
    stop(deparse("Select a data name"))
  
  
  for( i in 1:n.max){
    
    #Data Treatment
    if( i == 1){
      data <- as.matrix(data)
      k <- rep(NA,2) #only for don't stop the loop on the first time
    }
    if( i > 1 )
      data <- new.data
    pedigree <- asciitonumber(data,unk=unk)
    ind.data <- pedigree$ind.data
    sire <- pedigree$sire
    dire <- pedigree$dire
    ind <- c(1:length(sire))
    
    right.pos <- rep(NA, length=length(sire))
    
    #Verify alternatively sire and dire each iteraction+1
    parent <- sire
    parent.ind <- "sire"
    
    if( indicator%%2 == 1 ){
      parent <- dire
      parent.ind <- "dire"
    }
    
    
    for ( j in 1:length(parent)){
      if( is.na(match(parent[j], ind[1:j])) && parent[j] != 0){
        right.pos[j] <- which(ind == parent[j])
      }
    }
    
    error <- c()
    for( j in 1:length(parent))
      if( parent[j] > j)
        error <- c(error,j)
    
    #Print the step point
    if( save ){
      cat( paste("iteraction #",i,parent.ind,"\n",sep=""))
      cat( paste(error,"\n"))
    }
    
    
    #Right positions of the new data
    if( length(error) > 0 ){
      after <- which( !is.na(right.pos))
      before <- right.pos[after]
      for( j in 1:length(before)){
        if( j == 1 || before[j] != before[j-1] ){
          ind[after[j]] <- before[j]
          ind[before[j]] <- after[j]
        }
      }
    }
    new.data <- data[ind,]
    
    
    # Verify changes in the loop
    lastk <- k #indicator that there is no more changes in the last parent
    k <- length(error) #indicator that no changes
    
    
    if (  k == 0 )
      indicator = indicator+1
    
    if ( i != 1 && k == 0 && lastk == 0 ){
      cat("Your data was chronologically organized with success. \n")
      if(save){
        cat(paste("orgnew.txt",sep=""))
        write.table(new.data,  file=paste("orgped.txt",sep=""), quote=FALSE, row.names=FALSE, col.names=FALSE)
      }
      return(pedigree)
    }
    
    if ( i == n.max ){
      cat("Your data was not chronologically organized with sucess. Check your data with missing.data function and/or verify the individuals in the 2 last iteractions above descripted (the number is the row in the file: \n")
      cat(paste("orgped.txt",sep=""))
      write.table(new.data, file=paste("orgped.txt",sep=""), quote=FALSE, row.names=FALSE, col.names=FALSE)
      return(pedigree)
    }
  }
}

asciitonumber_legacy <- function(
    pedigree.data,
    unk=0
){
  if( ncol(pedigree.data) != 3 ){
    print("Data with more than 3 columns, please verify")
    return()
  }
  ind.data <- as.vector(c(unk,pedigree.data[,1]))
  sire.data <- as.vector(pedigree.data[,2])
  dire.data <- as.vector(pedigree.data[,3])
  sire <- match(sire.data, ind.data)
  dire <- match(dire.data, ind.data)
  ind <- as.vector(c(1:length(ind.data)))
  sire <- sire-1
  dire <- dire-1
  ind <- ind[-length(ind)]
  ind.data <- ind.data[-1]
  pedigree <- list(sire=sire,dire=dire,ind.data=pedigree.data[ind,1])
  return(pedigree)
}


sortped<-function(data = NULL, loop.in = 1000, loop.between = 100, print = FALSE)
{   
  if (is.null(data)) 
    stop(deparse("Please define the variable data"))
  
  stop.loop.1 <- stop.loop.2 <- FALSE
  
  for(j in 1:loop.between){
    if(print) cat(paste("looping between...",print(j)))
    
    for(i in 1:loop.in){
      if(print) cat(paste("looping in first parent...",print(i)))
      ind<-data[,1]
      sire<-data[,2]
      dire<-data[,3]
      index<-1:length(ind)
      compare<-match(dire,ind)
      compare[which(is.na(compare))]<-0
      loop<-which(compare>index)
      newindex<-index
      newindex[loop[1]]<-compare[loop[1]]
      newindex[compare[loop[1]]]<-loop[1]
      data<-data[newindex,]
      if(print) print(length(loop))
      if( length(loop) == 0 && i == 1)
        stop.loop.1 <- TRUE
      if(length(loop)==0) break
    }
    
    for(i in 1:loop.in){
      if(print) cat(paste("looping in second parent...",print(i))) 
      ind<-data[,1]
      sire<-data[,2]
      dire<-data[,3]
      index<-1:length(ind)
      compare<-match(sire,ind)
      compare[which(is.na(compare))]<-0
      loop<-which(compare>index)
      newindex<-index
      newindex[loop[1]]<-compare[loop[1]]
      newindex[compare[loop[1]]]<-loop[1]
      data<-data[newindex,]
      if(print) print(length(loop))
      if( length(loop) == 0 && i == 1)
        stop.loop.2 <- TRUE
      if( length(loop) == 0 ) break
    }
    
    if( stop.loop.1 && stop.loop.2) break
  }
  return(data)
}


test_that("Both versions of asciitonumber and datatreat give same output", {
  # Test data
  pedigree <- matrix(c(
    "A", "B", "C",
    "B", "D", "E",
    "C", "F", "G",
    "D", "0", "0",
    "E", "0", "0",
    "F", "0", "0",
    "G", "0", "0"
  ), ncol = 3, byrow = TRUE)
  
  # Set colnames (optional, for clarity)
  colnames(pedigree) <- c("id", "sire", "dam")
  
  # Run legacy R version
  ped_sorted_legacy <- sortped(pedigree)
  res_legacy <- asciitonumber_legacy(ped_sorted_legacy)
  
  # Run C++ version
  ped_sorted_cpp <- datatreat_cpp(pedigree, save = FALSE)
  res_cpp <- ascii_to_number(ped_sorted_legacy)
  
  # Compare sire
  expect_equal(res_cpp$sire, res_legacy$sire)
  
  # Compare dam
  expect_equal(res_cpp$dire, res_legacy$dire)
  
  # Compare individual names
  expect_equal(as.character(res_cpp$ind_data), as.character(res_legacy$ind.data))
})

test_that("Second test for pedigree", {
  pedigree <- matrix(c(
    "A", "B", "C",
    "B", "D", "E",
    "C", "F", "G",
    "D", "0", "0",
    "E", "0", "0",
    "F", "0", "0",
    "G", "0", "0"
  ), ncol = 3, byrow = TRUE)
  colnames(pedigree) <- c("id", "sire", "dam")
  
  # Legacy
  ped_sorted_legacy <- sortped(pedigree)
  res_legacy <- asciitonumber_legacy(ped_sorted_legacy)
  
  # C++
  ped_sorted_cpp <- datatreat_cpp(pedigree, save = FALSE)
  res_cpp <- ascii_to_number(ped_sorted_legacy)
  
  expect_equal(res_cpp$sire, res_legacy$sire)
  expect_equal(res_cpp$dire, res_legacy$dire)
  expect_equal(as.character(res_cpp$ind_data), as.character(res_legacy$ind.data))
})


test_that("Handles missing parents as '0' or NA", {
  pedigree <- matrix(c(
    "A", "B", "C",
    "B", "D", NA,
    "C", "0", "G",
    "D", NA, "0",
    "E", "0", "0",
    "F", NA, NA,
    "G", "0", "0"
  ), ncol = 3, byrow = TRUE)
  colnames(pedigree) <- c("id", "sire", "dam")
  
  # Convert NAs to "0" to match unknown handling
  pedigree[is.na(pedigree)] <- "0"
  
  ped_sorted_legacy <- sortped(pedigree)
  res_legacy <- asciitonumber_legacy(ped_sorted_legacy)
  ped_sorted_cpp <- datatreat_cpp(pedigree, save = FALSE)
  res_cpp <- ascii_to_number(ped_sorted_legacy)
  
  expect_equal(res_cpp$sire, res_legacy$sire)
  expect_equal(res_cpp$dire, res_legacy$dire)
  expect_equal(as.character(res_cpp$ind_data), as.character(res_legacy$ind.data))
})


test_that("Unordered pedigree gets sorted correctly", {
  pedigree <- matrix(c(
    "E", "F", "G",
    "A", "B", "C",
    "B", "0", "0",
    "C", "0", "0",
    "F", "0", "0",
    "D", "E", "0",
    "G", "0", "0"
  ), ncol = 3, byrow = TRUE)
  colnames(pedigree) <- c("id", "sire", "dam")
  
  ped_sorted_legacy <- sortped(pedigree)
  res_legacy <- asciitonumber_legacy(ped_sorted_legacy)
  ped_sorted_cpp <- datatreat_cpp(pedigree, save = FALSE)
  res_cpp <- ascii_to_number(ped_sorted_legacy)
  
  expect_equal(res_cpp$sire, res_legacy$sire)
  expect_equal(res_cpp$dire, res_legacy$dire)
  expect_equal(as.character(res_cpp$ind_data), as.character(res_legacy$ind.data))
})


test_that("Large pedigree works the same", {
  set.seed(42)
  ids <- paste0("ID", 1:100)
  parents <- sample(c(ids, "0"), size = 200, replace = TRUE)
  pedigree <- matrix(cbind(ids, parents[1:100], parents[101:200]), ncol = 3)
  colnames(pedigree) <- c("id", "sire", "dam")
  
  # Clean unknowns if needed
  pedigree[is.na(pedigree)] <- "0"
  
  ped_sorted_legacy <- sortped(pedigree)
  res_legacy <- asciitonumber_legacy(ped_sorted_legacy)
  ped_sorted_cpp <- datatreat_cpp(pedigree, save = FALSE)
  res_cpp <- ascii_to_number(ped_sorted_legacy)
  
  expect_equal(res_cpp$sire, res_legacy$sire)
  expect_equal(res_cpp$dire, res_legacy$dire)
  expect_equal(as.character(res_cpp$ind_data), as.character(res_legacy$ind.data))
})


test_that("Missing parent not in ID column still handled", {
  pedigree <- matrix(c(
    "A", "X", "Y",
    "B", "A", "Z",
    "C", "B", "0"
  ), ncol = 3, byrow = TRUE)
  colnames(pedigree) <- c("id", "sire", "dam")
  pedigree[is.na(pedigree)] <- "0"
  
  ped_sorted_legacy <- sortped(pedigree)
  res_legacy <- asciitonumber_legacy(ped_sorted_legacy)
  ped_sorted_cpp <- datatreat_cpp(pedigree, save = FALSE)
  res_cpp <- ascii_to_number(ped_sorted_legacy)
  
  expect_equal(res_cpp$sire, res_legacy$sire)
  expect_equal(res_cpp$dire, res_legacy$dire)
  expect_equal(as.character(res_cpp$ind_data), as.character(res_legacy$ind.data))
})

test_that("Duplicated IDs are caught or warned in the new C++", {
  pedigree <- matrix(c(
    "A", "B", "C",
    "B", "D", "E",
    "A", "F", "G"  # duplicate "A"
  ), ncol = 3, byrow = TRUE)
  colnames(pedigree) <- c("id", "sire", "dam")
  
  # C++ version should stop with an error
  expect_error({
    datatreat_cpp(pedigree)
  }, regexp = "Duplicate individual ID found", ignore.case = TRUE)
})

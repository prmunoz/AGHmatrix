#########################################################################
# 									
# Package: AGHmatrix 							
# 									
# File: AmatrixPolyCross.R
# Contains: AmatrixPolyCross 							
# 									
# Written by Rodrigo Rampazo Amadeu 			
# Contributors: Leticia AC Lara
# 									
# First version: 20-Nov-2020 
# Last update: 05-Aug-2021 						
# License: GPL-3
# 								
#########################################################################

#' Construction of pedigree-based relationship matrix with parental guessing possibility
#'
#' Creates an additive relationship matrix A based on a non-deterministic pedigree with 4+ columns where each column represents a possible parent. This function was built with the following designs in mind. 
#' 1) A mating design where you have equally possible parents. For example, a generation of insects derived from the mating of three insects in a cage. All the insects in this generation will have the same expected relatedness with all the possible parents (1/3). If there are only two parents in the cage, the function assumes no-inbreeding and the pedigree is deterministic (the individual is offspring of the cross between the two parents). Another example, a population of 10 open-pollinated plants where you harvest the seeds without tracking the mother. 
#' 2) When fixedParent is TRUE: a mating design where you know one parent and might know the other possible parents. For example, a polycross design where you have seeds harvested from a mother plant and possible polen donors.
#'
#' @param data pedigree data name. Unknown value should be equal 0. See example for construction.
#' @param fixedParent if false, assumes that all the parents are equally possible parents. If true, assumes that the first parental is known and the others are equally possible parents. Default = FALSE.
#' 
#' @return Matrix with the relationship between the individuals.
#'
#' @examples
#' #the following pedigree has the id of the individual followed by possible parents
#' #if 0 is unknown
#' #the possible parents are filled from left to right
#' #in the pedigree data frame examples:
#' #id 1,2,3,4 have unknown parents and are assumed unrelated
#' #id 5 has three possible parents (1,2,3)
#' #id 6 has three possible parents (2,3,4)
#' #id 7 has two parents (deterministic case here, the parents are 3 and 4)
#' #id 8 has four possible parents (5,6,7,1)
#' 
#' pedigree = data.frame(id=1:8,
#'                       parent1 = c(0,0,0,0,1,2,3,5),
#'                       parent2 = c(0,0,0,0,2,3,4,6),
#'                       parent3 = c(0,0,0,0,3,4,0,7),
#'                       parent4 = c(0,0,0,0,0,0,0,1),
#'                       parent5 = 0)
#'
#' print(pedigree)
#'
#' AmatrixPolyCross(pedigree)
#' 
#' #when polyCross is set to be true:
#' #id 5 is offspring of parent 1 in a deterministic way and two other possible parents (2,3)
#' #id 6 is offspring of parent 2 in a deterministic way and two other possible parents (3,4)
#' #id 7 has two parents (deterministic case here, the parents are 3 and 4); as before
#' #id 8 is offspring of parent 5 in a deterministic way and has three other possible parents (6,7,1)
#' 
#' AmatrixPolyCross(pedigree,fixedParent=TRUE)
#'
#' @author Rodrigo R Amadeu, \email{rramadeu@@gmail.com}
#'
#' @export
#' 
AmatrixPolyCross = function(data = NULL, fixedParent=FALSE){
  unk = 0
  orig.order <- as.character(data[, 1])
  
  data1 <- AGHmatrix::datatreat(data = data[,c(1,2,3)], unk = unk)
  Parents <- data1$sire
  Parents <- rbind(Parents,data1$dire)
  
  for(i in 4:ncol(data))
    Parents <- rbind(Parents, AGHmatrix::datatreat(data = data[,c(1,2,i)], unk = unk)$dire)
  
  for(i in 1:ncol(Parents)){
    tst = which(Parents[,i]!=0)
    if(length(tst)>0)
      if(min(which(Parents[,i]==0)) < max(tst))
        stop(deparse(paste0("Check parent order of line ",i,
                            ", missing (or non-used) values should be located at the right!")))
  }
  
  if (length(data$sire) > 1000) 
    cat("Processing a large pedigree data... It may take a couple of minutes... \n")
  
  if(fixedParent==FALSE){
    n <- ncol(Parents)
    Time = proc.time()
    
    combs = NULL
    for(i in 2:nrow(Parents)){
      combs = cbind(combs, combn(i,2))
    }
    
    combs = t(unique(t(combs))) #fixing order by number of max parents
    
    cat("Constructing matrix A using ploidy = 2 \n")
    Afinal = matrix(NA,n,n)
    Afinal[1, 1] <- 1
    Acombs <- array(NA, dim=c(n,n,ncol(combs)))
    
    for (i in 2:n) {
      for(combsIndex in 1:ncol(combs)){
        A <- Afinal
        s <- Parents[combs[1,combsIndex],]
        d <- Parents[combs[2,combsIndex],]
        if (s[i] == 0 && d[i] == 0) {
          A[i, i] <- 1
          for (j in 1:(i - 1)) A[j, i] <- A[i, j] <- 0
        }
        if (s[i] == 0 && d[i] != 0) {
          A[i, i] <- 1
          for (j in 1:(i - 1)) A[j, i] <- A[i, j] <- 0.5 * 
              (A[j, d[i]])
        }
        if (d[i] == 0 && s[i] != 0) {
          A[i, i] <- 1
          for (j in 1:(i - 1)) A[j, i] <- A[i, j] <- 0.5 * 
              (A[j, s[i]])
        }
        if (d[i] != 0 && s[i] != 0) {
          A[i, i] <- 1 + 0.5 * (A[d[i], s[i]])
          for (j in 1:(i - 1)) A[j, i] <- A[i, j] <- 0.5 * 
              (A[j, s[i]] + A[j, d[i]])
        }
        Acombs[1:i,1:i,combsIndex] = A[1:i,1:i]
      }
      
      ## Total parents
      totalParents = length(which(Parents[,i]!=0))
      if(totalParents < 3){
        Afinal[i,] = Afinal[,i] = Acombs[i,,1]
      }else{
        Afinal[i,] = Afinal[,i] = apply(Acombs[i,,(1:choose(totalParents,2))],1,mean)
      }
    }
    A = Afinal
  }
  
  if(fixedParent==TRUE){
    n <- ncol(Parents)
    Time = proc.time()
    Afinal = matrix(NA,n,n)
    Afinal[1, 1] <- 1
    Acombs <- array(NA, dim=c(n,n,(nrow(Parents)-1)))
    s <- Parents[1,]   
    for (i in 2:n) {
      for(combsIndex in 1:(nrow(Parents)-1)){
        A <- Afinal
        d <- Parents[combsIndex+1,]
        if (s[i] == 0 && d[i] == 0) {
          A[i, i] <- 1
          for (j in 1:(i - 1)) A[j, i] <- A[i, j] <- 0
        }
        if (s[i] == 0 && d[i] != 0) {
          A[i, i] <- 1
          for (j in 1:(i - 1)) A[j, i] <- A[i, j] <- 0.5 * 
              (A[j, d[i]])
        }
        if (d[i] == 0 && s[i] != 0) {
          A[i, i] <- 1
          for (j in 1:(i - 1)) A[j, i] <- A[i, j] <- 0.5 * 
              (A[j, s[i]])
        }
        if (d[i] != 0 && s[i] != 0) {
          A[i, i] <- 1 + 0.5 * (A[d[i], s[i]])
          for (j in 1:(i - 1)) A[j, i] <- A[i, j] <- 0.5 * 
              (A[j, s[i]] + A[j, d[i]])
        }
        Acombs[1:i,1:i,combsIndex] = A[1:i,1:i]
      }
      
      ## Total parents
      totalParents = length(which(Parents[,i]!=0))
      if(totalParents < 3){
        Afinal[i,] = Afinal[,i] = Acombs[i,,1]
      }else{
        Afinal[i,] = Afinal[,i] = apply(Acombs[i,,1:(totalParents-1)],1,mean)
      }
    }
    A = Afinal
  }
  
  Time = as.matrix(proc.time() - Time)
  cat("Completed! Time =", Time[3]/60, " minutes \n")
  rownames(A) <- colnames(A) <- data1$ind.data
  A <- A[orig.order, orig.order]
  return(A)
}



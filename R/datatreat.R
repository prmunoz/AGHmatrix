#########################################
# 									
# Package: AGHmatrix 							
# 									
# File: datatreat.R							
# Contains: datatreat asciitonumber 
# 									
# Written by Rodrigo Rampazo Amadeu 					
# 									
# First version: Feb-2014 					
# Last update: 24-Apr-2015 						
# License: GPL-3	
# 									
#########################################

#' Organizes pedigree data in a chronological way
#'
#' This function organizes pedigree data in a chronological way and return 3 lists: i) parental 1 values (numeric); ii) parental 2 values (numeric); iii) real names of the individuals. Also save a .txt file with new pedigree file.
#' @param data name of the pedigree data frame. Default=NULL.
#' @param unk the code of the data missing. Default=0.
#' @param n.max max number of iteractions to get the chronological order. Default = 50
#' @param save if TRUE, save the genealogy in a .txt file
#'
#' @return list with parental 1, parental 2, and real names of the individuals (key) also saves a txt file with the new chronological pedigree.
#'
#' @examples 
#' data(ped.mrode)
#' datatreat(ped.mrode)
#'
#' @author Rodrigo R Amadeu, \email{rramadeu@@gmail.com}
#'
#' @export

datatreat <- function(data=NULL,
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

# This function creates a list with numeric indices given a pedigree data. 
# Also it checks if all the listed parent name are listed before in the individual 
# name column and if the parent exist in the matrix . 
asciitonumber <- function(
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

# This function organizes pedigree data in a chronological way and return 3 lists: 
# i) parental 1 values (numeric); ii) parental 2 values (numeric); iii) real names of 
# the individuals. Also save a .txt file with new pedigree file.
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

# This function verify which rows in a pedigree data has missing parental or conflictuos data
verifyped <- function(pedigree,
                      unk=0
){
  flag<-FALSE
  if( ncol(pedigree) != 3 ){
    print("Data with more than 3 columns, please verify")
    flag<-TRUE
    return(flag)
  }
  
  if(length(unique(pedigree[,1]))<nrow(pedigree)){
    print("Data with repeated entry, please verify the following entries lines")
    print(which(duplicated(pedigree[,1]),arr.ind=TRUE))
    flag<-TRUE
    return(flag)
  }
  
  #Treating all as numeric
  ind.data <- as.vector(pedigree[,1])
  sire.data <- as.vector(pedigree[,2])
  dire.data <- as.vector(pedigree[,3])
  sire <- match(sire.data, c(ind.data,"0"))
  dire <- match(dire.data, c(ind.data,"0"))
  ind <- as.vector(c(1:length(ind.data)))
  
  missing <- c()   
  #Verify the individual w/ same name in sire/dire
  missing$conflict <- c(which(sire == ind),which( dire == ind))
  if(length(missing$conflict)>0){
    print("The following rows have the individual name equals to the parental name. Please verify.")
    print(pedigree[missing$conflict,])
    flag<-TRUE
  }
  
  #Verify the missing sire (Parent 1)
  missing$sire.na <- c( which(is.na(sire)))
  if(length(missing$sire)>0){
    print("The following rows have the parental 1 name (column 2) missing in the pedigree. Please verify.")
    print(pedigree[missing$sire,])
    flag<-TRUE
  }
  
  #Verify the missing dire (Parent 2)
  missing$dire.na <- c( which(is.na(dire)))
  if(length(missing$dire)>0){
    print("The following rows have the parental 2 name (column 3) missing in the pedigree. Please verify.")
    print(pedigree[missing$dire,])
    flag<-TRUE
  }
  return(flag)
}

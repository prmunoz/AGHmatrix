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

#' Organizes pedigree data in a chronological way using C++
#'
#' This function organizes pedigree data in a chronological way and returns a list:
#' i) parental 1 values (numeric); ii) parental 2 values (numeric); iii) individual names.
#' Also saves a `.txt` file if `save = TRUE`.
#'
#' @param data A 3-column data frame or matrix with individual, sire, and dam
#' @param unk Code for unknown parent (default = 0)
#' @param n.max Max number of iterations (default = 50)
#' @param save If TRUE, saves the reordered pedigree to a text file
#'
#' @return A list with elements: sire (numeric), dire (numeric), ind_data (names)
#' @author Rodrigo R. Amadeu \email{rramadeu@gmail.com}
#' @author Thiago de Paula Oliveira
#' 
#' @export
#'
#' @examples
#' # data(ped.mrode)
#' # datatreat(ped.mrode)
datatreat <- function(data = NULL,
                      n.max = 50,
                      unk = 0,
                      save = FALSE) {
  if (is.null(data)) {
    stop("Select a data name")
  }
  data <- as.matrix(data)
  datatreat_cpp(data, n_max = n.max, unk = as.character(unk), save = save)
}

# This function creates a list with numeric indices given a pedigree data. 
# Also it checks if all the listed parent name are listed before in the individual 
# name column and if the parent exist in the matrix . 
#' Convert pedigree names to numeric indices using C++
#'
#' @param pedigree.data A 3-column matrix or data frame with individual, sire, and dam
#' @param unk Code for unknown parent (default = "0")
#' @return A list with sire, dam, and individual names
#' @author Rodrigo R. Amadeu \email{rramadeu@gmail.com}
#' @author Thiago de Paula Oliveira
#' @export
asciitonumber <- function(pedigree.data, unk = 0) {
  pedigree.data <- as.matrix(pedigree.data)
  if (ncol(pedigree.data) != 3) {
    stop("Data with more than 3 columns, please verify")
  }
  ascii_to_number(pedigree.data, as.character(unk))
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

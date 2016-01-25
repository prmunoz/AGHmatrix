#########################################################################
# 									#
# Package: AGHmatrix 							#
# 									#
# File: verifyped.R 							#
# Contains: verifyped 							#
# 									#
# Written by Rodrigo Rampazo Amadeu 					#
# 									#
# First version: Feb-2014 						#
# Last update: 15-Dec-2015 						#
# License: GNU General Public License version 2 (June, 1991) or later 	#
# 									#
#########################################################################

#' Survying on verify pedigree
#'
#' This function verify which rows in a pedigree data has missing parental or conflictuos data
#'
#' @param pedigree data name from a pedigree list
#' @param unk unknown value of your data
#'
#' @return i) warnings with the conflicting data regarding number of columns in the data, repeated entries, missing parental names in the entries. ii) TRUE if there was any warning.
#'
#' @examples verifyped()
#'
#' @author Rodrigo R Amadeu, \email{rramadeu@@gmail.com}
#'
#' @export


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


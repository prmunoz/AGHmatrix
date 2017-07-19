#########################################
# 									
# Package: AGHmatrix 							
# 									
# File: missingdata.R							
# Contains: missingdata							
# 									
# Written by Rodrigo Rampazo Amadeu 					
# 									
# First version: Feb-2014 						
# Last update: 14-Apr-2015 						
# License: GPL-3
# 									
#########################################

#' Survying on missing data
#'
#' This function verify which rows in a pedigree data has missing parental or conflictuos data
#'
#' @param data data name from a pedigree list
#' @param unk unknown value of your data
#'
#' @return list with $conflict: rows of the data which are at least one parental name equal to the individual. $missing.sire: rows of the data which arie missing data sire (Parental 1) information. $missing.dire: same as above for dire (Parental 2). $summary.missing: summary of the missing data. 2 columns, 1st for the name of the parental listed, 2nd for the how many times appeared in the data.
#'
#' @examples 
#' data(ped.mrode)
#' missingdata(ped.mrode)
#'
#' @author Rodrigo R Amadeu, \email{rramadeu@@gmail.com}
#'
#' @export

missingdata <- function(data,
			unk=0
                        ){
    pedigree.data<-data
    data <- c()

    if( ncol(pedigree.data) != 3 ){
      print("Data with more than 3 columns, please verify")
      return()
    }

    #Treating all as numeric
    ind.data <- as.vector(c(unk,as.character(pedigree.data[,1])))
    sire.data <- as.vector(pedigree.data[,2])
    dire.data <- as.vector(pedigree.data[,3])
    sire <- match(sire.data, ind.data)
    dire <- match(dire.data, ind.data)
    ind <- as.vector(c(1:length(ind.data)))
    missing <- c()

    #Verify the individual w/ same name in sire/dire
    missing$conflict <- c(which( sire == ind[-1]),which( dire == ind[-1] ))

    #Verify the missing sire (Parent 1)
    missing$sire.na <- c( which(is.na(sire)))

    #Verify the missing dire (Parent 2)
    missing$dire.na <- c( which(is.na(dire)))

    #Making a summary of the missing data
    missing$sire <- as.matrix(summary(as.factor(pedigree.data[which(is.na(sire)),2])))
    missing$dire <- as.matrix(summary(as.factor(pedigree.data[which(is.na(dire)),3])))
    names <- unique(c(rownames(missing$sire),rownames(missing$dire)))
    pos.sire <- match(rownames(missing$sire),names)
    pos.dire <- match(rownames(missing$dire),names)
    missing$parent <- rep(0,length(names))
    missing$parent[pos.dire] <- missing$dire
    missing$parent[pos.sire] <- missing$parent[pos.sire]+missing$sire
    missing$parent <- as.matrix(missing$parent)
    rownames(missing$parent) <- names
    #Final list
    missing <- list(conflict=missing$conflict,missing.sire=missing$sire.na,missing.dire=missing$dire.na,summary.missing=missing$parent)

    ## Improve this part :)
    #if(molecular){
    # if( csv )
    #    data <- read.csv("molecular_diploid.csv")
    #  mol.data <- data[,-1]
    #  row.names(mol.data) <- data[,1]
    #  mol.data <- replace(mol.data, mol.data == unk, NA)

    #  missing.per.ind <- apply(mol.data,1, function(x) sum(is.na(x)))# / ncol(example) * 100
    #  missing.per.marker  <- apply(t(mol.data),1,function(x) sum(is.na(x)))
    #  summary(missing.per.ind)/length(missing.per.marker)
    #  summary(missing.per.marker)/length(missing.per.ind)
    #}
  return(missing)
}


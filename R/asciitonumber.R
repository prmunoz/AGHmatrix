#########################################################################
# 									#
# Package: AGHmatrix 							#
# 									#
# File: asciitonumber.R							#
# Contains: asciitonumber						#
# 									#
# Written by Rodrigo Rampazo Amadeu 					#
# 									#
# First version: Feb-2014 						#
# Last update: 14-Apr-2015 						#
# License: GNU General Public License version 2 (June, 1991) or later 	#
# 									#
#########################################################################

#' Given a pedigree ascii format returns a numeric list
#'
#' This function creates a list with numeric indices given a pedigree data. Also it checks if all the listed parent name are listed before in the individual name column and if the parent exist in the matrix . 
#'
#' @param pedigree.data matrix with 3 columns (without row names/column names), the first for the individual name, second for parent 1, third for parent 2.
#' @param unk code of the data missing, default = 0
#' @param force replace all the missing parents with the unk value
#'
#' @return List with $sire: vector with the numbers respectives to the first column parent in the data; $dire: as sire but for second column; $ind.data: the matrix of original data.
#' 
#' @examples asciitonumber()
#'
#' @author Rodrigo R Amadeu, \email{rramadeu@@gmail.com}
#'
#' @export


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

#########################################################################
# 									
# Package: AGHmatrix 							
# 									
# File: expandAmatrix.R
# Contains: expandAmatrix 							
# 									
# Written by Rodrigo Rampazo Amadeu
# 									
# First version: Oct-2021 
# Last update: 03-Nov-2019 						
# License: GPL-3
# 								
#########################################################################

#' Add new crosses to a current A matrix
#'
#' Expand a current A matrix with a new pedigree. The parents in the new pedigree should also be in the A matrix.
#'
#' @param newPedigree pedigree data name (3-column way format). Unknown value should be equal 0.
#' @param A numerator relationship matrix output from Amatrix function.
#' @param returnAll if TRUE returns old A with new A, if FALSE returns only new A
#' 
#' @return Matrix with the Relationship between the individuals.
#'
#' @examples
#' data(ped.sol)
#' ped.initial = ped.sol[1:1120,]
#' ped.new = ped.sol[-c(1:1120),]
#' #Computing additive relationship matrix:
#' A = Amatrix(ped.initial, ploidy=2)
#' Anew = expandAmatrix(ped.new, A)
#' 
#' #Comparing with one-step building..
#' Afull = Amatrix(ped.sol, ploidy=2)
#' test = Anew-Afull
#' which(test!=0)
#' 
#' @author Rodrigo R Amadeu, \email{rramadeu@@gmail.com}
#' 
#' @export

expandAmatrix <- function(newPedigree = NULL,
                          A=NULL,
                          returnAll=TRUE){

  if( is.null(newPedigree))
    stop(deparse("newPedigree argument is missing"))
  
  if( is.null(A))
    stop(deparse("A argument is missing"))
  
    
  if(any(is.na(match(newPedigree[,2],rownames(A))))+any(is.na(match(newPedigree[,3],rownames(A)))))
    stop(deparse("There are individuals in the new pedigree with missing parents in A"))
  
  ## Creating a line of 0s
  A = rbind(A,rep(0,nrow(A)))
  A = cbind(A,c(rep(0,ncol(A)),1))
  colnames(A)[ncol(A)] = rownames(A)[ncol(A)] = "0"
  
  NewA = A[match(newPedigree[,2],rownames(A)),match(newPedigree[,2],rownames(A))] +
    A[match(newPedigree[,2],rownames(A)),match(newPedigree[,3],rownames(A))] + 
    A[match(newPedigree[,3],rownames(A)),match(newPedigree[,2],rownames(A))] + 
    A[match(newPedigree[,3],rownames(A)),match(newPedigree[,3],rownames(A))]
  
  NewA = NewA/4
  diag(NewA) = 1+diag(A[match(newPedigree[,2],rownames(A)),match(newPedigree[,3],rownames(A))])/2
  
  if(!returnAll){
    rownames(NewA)=colnames(NewA)=newPedigree[,1]
    return(NewA)
  }else{
    NewOldA = A[match(newPedigree[,2],rownames(A)),] + A[match(newPedigree[,3],rownames(A)),]
    NewOldA = NewOldA/2
    NewOldA = NewOldA[,-ncol(NewOldA)]
    NewOldA = cbind(rbind(A[-ncol(A),-ncol(A)],NewOldA),
                    rbind(t(NewOldA),NewA))
    return(NewOldA)
  }
}

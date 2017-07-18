#########################################################################
# 									
# Package: AGHmatrix 							
# 									
# File: formatmatrix.R 						
# Contains: formatmatrix						
# 									
# Written by Rodrigo Rampazo Amadeu 					
# 									
# First version: Feb-2014 						
# Last update: 20-May-2015 						
# License: GPL-3
# 								
#######################################

#' Transform a matrix in 3 columns
#'
#' Given any square matrix transform it in a 3 columns way (row, column, value) mainly to be used in outsourcing data processing (as ASREML-standalone)
#'
#' @param data matrix (nxn).
#' @param save if TRUE save the output in a file. Default=TRUE.
#' @param return if TRUE return the output in a object. Default=FALSE.
#' @param name name of the csv file to be saved. Default=data name.
#' @param round.by select the number of digits after 0 you want in your data. Default = 12
#' @param exclude.0 if TRUE, remove all lines equal to zero (ASREML option). Default = TRUE
#' 
#' @return a object or a csv file with a table with 3 columns representing the matrix.
#'
#' @examples
#' #Example with random matrix
#' data<-matrix(c(1,0.1,0,0.1,1,0,0,0,1.1),3)
#' formatmatrix(data=data,save=FALSE,return=TRUE,exclude.0=TRUE)
#'
#' #Example with pedigree matrix
#' #Reading the example data
#' data(ped.mrode)
#' #Making Relationship Matrix
#' Amrode<-Amatrix(ped.mrode)
#' #Inverting the Matrix
#' Amrode.inv<-solve(Amrode)
#' #Making the 3 columns format
#' Amrode.inv.ASREML<-formatmatrix(Amrode,save=FALSE,return=TRUE,exclude.0=TRUE)
#' #Printing it
#' Amrode.inv.ASREML 
#'
#' @author Rodrigo R Amadeu, \email{rramadeu@@gmail.com}
#'
#' @export

formatmatrix <- function(
    data = NULL,
    save = TRUE,
    return = FALSE,
    name = deparse(substitute(data)),
    round.by = 12,
    exclude.0 = TRUE)
    {
  if( is.null(data))
    stop(deparse("Define a matrix data"))
  
  cat("Converting to column format... \n")
  
  Time <- proc.time()
  n <- nrow(data)
  first <- second <- third <- c()
  for( i in 1:n ){
    first <- c(first,rep(i,i))
    second <- c(second,c(1:i))
  }
  third <- data[upper.tri(data,diag=TRUE)]
  
  third<-round(third,round.by)
  
  columns <- cbind(first, second, third)
  
  if(exclude.0)
    columns<-columns[third!=0,]
  
  if(save){
    write.table(columns,file=paste(name,".csv",sep=""),sep=" ",row.names=FALSE,col.names=FALSE)
    cat(paste("Saved as ",name,".csv"," \n",sep=""))
  }
  
  Time = as.matrix(proc.time()-Time[3])
  
  cat("Completed! Time =", Time[3]/60," minutes \n")
  if(return)
    return( columns )
}

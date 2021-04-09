#########################################
# 									
# Package: AGHmatrix 							
# 									
# File: filterpedigree.R							
# Contains: filterpedigree
# 									
# Written by Rodrigo Rampazo Amadeu 					
# 									
# First version: Apr-2021 					
# Last update: 09-Apr-2021 						
# License: GPL-3	
# 									
#########################################

#' Filter the pedigree to keep only the genealogy of a subset of individuals
#'
#' Filter the pedigree to keep only the genealogy of a subset of individuals
#' @param inds vector with strings of individuals to keep their genealogy in the matrix
#' @param data name of the pedigree data frame. Default=NULL.
#'
#' @return a data frame with pedigree containing the genealogy of the selected individuals
#' 
#' @examples 
#' data(ped.sol)
#' new.ped.sol = filterpedigree(inds = c("MSW168-2","W14090-3","W14090-4"),data=ped.sol)
#'
#' @author Rodrigo R Amadeu, \email{rramadeu@@gmail.com}
#'
#' @export

filterpedigree <- function(inds, data){
  output <- NULL
  progress <- round(length(inds)/10)
  perc <- 10
  for(i in 1:length(inds)){
    if(length(inds)>100){
      if(i %% progress ==0){
        cat(paste0(perc,"% \n"))
        perc=perc+10
      }
    }
    ped_out<- data[which(data[,1] == inds[i]),]
    if(nrow(ped_out)==0){
      stop(deparse(paste(inds[i],"doesn't exist in this pedigree.")))
    }
    trigger <- 1
    while(trigger>0){
      init <- nrow(ped_out)
      ped_in <- data[which(data[,1] %in% c(ped_out[,2],ped_out[,3])),]
      ped_out <- unique(rbind(ped_in, ped_out))
      trigger <- nrow(ped_out)-init
    }
    output = rbind(output,ped_out)
    output = unique(output)
  }
  return(output)
}
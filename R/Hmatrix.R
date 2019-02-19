########################################
# 									
# Package: AGHmatrix 							
# 									
# File: Hmatrix.R 							
# Contains: Hmatrix 							
# 									
# Written by Rodrigo Rampazo Amadeu 					
# 									
# First version: Feb-2014 						
# Last update: 22-Apr-2015 						
# License: GPL-3
# 									
#######################################

#' Construction of Relationship Matrix H
#'
#' Given a matrix A and a matrix G returns a H matrix - A matrix corrected by G.
#'
#' @param A A matrix from function Amatrix
#' @param G G matrix from function Gmatrix
#' @param markers matrix marker which generated Gmatrix
#' @param c constant value of H computation, default: c=0
#' @param explore if TRUE performs exploratory analysis of the matrix
#' @param missingValue missing value in data. Default=-9.
#' @param maf max of missing data accepted to each marker. Default=0.05.
#' @param ploidy data ploidy (an even number between 2 and 20). Default=2.
#'
#' @return H Matrix with the relationship between the individuals based on pedigree and corrected by molecular information
#'
#' @examples 
#' data(ped.mrode)
#' #Build Amatrix diploid (no double reduction proportion)
#' Amat <- Amatrix(data=ped.mrode,ploidy=2)
#' markers <- matrix(c(0,0,0,0, 2,2,1,1, 1,1,0,1, 1,1,2,0, 2,1,1,0, 2,0,1,2),nrow=6, byrow=TRUE)
#' rownames(markers) <- rownames(Amat)
#' Gmat <- Gmatrix(markers)
#' Hmatrix(Amat,Gmat,markers)
#'
#' @author Rodrigo R Amadeu, \email{rramadeu@@gmail.com}
#' 
#' @references \emph{Munoz, P. R., Resende, M. F. R., Gezan, S. A., Resende, M. D. V., de los Campos, G., Kirst, M., Huber, D., Peter, G. F. (2014). Unraveling additive from nonadditive effects using genomic relationship matrices. Genetics, 198.4: 1759-1768.}
#'
#' @export

Hmatrix <- function(A=NULL,
                    G=NULL,
                    markers=NULL,
                    c=0,
                    explore=FALSE,
                    missingValue=-9,
                    maf=0,
                    ploidy=2
                    ){
    Aorig <- A
    Gorig <- G
    markersmatrix <- Gmatrix(markers,method="MarkersMatrix",ploidy=ploidy,missingValue=missingValue,maf=maf)

    Time = proc.time()
    cat("Comparing the matrices... \n")
    An <- rownames(Aorig)
    Gn <- rownames(Gorig)
    missingGmatrix <- which(is.na(match(An,Gn)))
    missingAmatrix <- which(is.na(match(Gn,An)))
    if(length(missingAmatrix)>0){
      Gnhat <- Gn[-missingAmatrix]
    }else{
      Gnhat <- Gn
    }
    
    if(length(missingGmatrix)>0){
      Anhat <- An[-missingGmatrix]    
    }else{
      Anhat <- An
    }
  
    A <- Aorig[Anhat,Anhat]
    G <- Gorig[Gnhat,Gnhat]

    missingGmatrix <- An[missingGmatrix]
    missingAmatrix <- Gn[missingAmatrix]

    Time = as.matrix(proc.time()-Time)
    cat("Completed! Time =", Time[3]/60," minutes \n")

    cat("Computing the H matrix... \n")

    Time = proc.time()
   #Computing the Variance of G by A classes
    classes <- as.numeric(levels(as.factor(A)))
    n <- length(classes)
    varA <- meanG <- matrix(NA,nrow=nrow(A),ncol=nrow(A))
    varAclasses <- c()
    for(i in 1:n){
        varA[A==classes[i]] <- varAclasses[i] <- var(G[A==classes[i]])
        meanG[A==classes[i]]  <- mean(G[A==classes[i]])
    }
    varA[is.na(varA)]<-mean(varAclasses,na.rm=TRUE) 

    #Computaing beta and H
    beta <- 1 - (c+(1/(markersmatrix[Gnhat,Gnhat]))/varA)
    H <- beta*(G-A)+A ######
    Aorig[Anhat,Anhat] = H
    cat("\n",length(missingGmatrix),"Individuals in A but not in G:",missingGmatrix,"\n")
    cat("\n",length(missingAmatrix),"Individuals in G but not in A:",missingAmatrix,"\n")
    Time = as.matrix(proc.time()-Time)
    cat("\n","Completed! Time =", Time[3]/60," minutes \n")
    cat("\n","Returning H = A matrix corrected by G... \n")
    return(Aorig)
}

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
# Last update: 09-Jul-2019 						
# License: GPL-3
# 									
#######################################

#' Construction of Combined Relationship Matrix H
#'
#' Given a matrix A and a matrix G returns a H matrix.
#'
#' @param A A matrix from function Amatrix
#' @param G G matrix from function Gmatrix
#' @param markers matrix marker which generated the Gmatrix
#' @param c constant value of H computation, default: c=0
#' @param method "Martini" or "Munoz", default="Martini"
#' @param missingValue missing value in data, default=-9.
#' @param maf max of missing data accepted to each markerm default=0.05.
#' @param ploidy data ploidy (an even number between 2 and 20), default=2.
#' @param tau to be used for Martini's method, default=1. 
#' @param omega to be used of Martini's method, default=1.
#' @param roundVar Munoz's method, how many digits to consider the relationship be of same class, default=2.
#' 
#' @return H Matrix with the relationship between the individuals based on pedigree and corrected by molecular information
#'
#' @examples 
#' \dontrun{
#' data(ped.sol)
#' data(snp.sol)
#' #Computing the numerator relationship matrix 10% of double-reduction
#' Amat <- Amatrix(ped.sol, ploidy=4, w = 0.1)
#' #Computing the additive relationship matrix based on VanRaden (modified)
#' Gmat <- Gmatrix(snp.sol, ploidy=4, missingValue=-9,
#'                 maf=0.05, method="VanRaden")
#' #Computing H matrix (Martini)
#' Hmat_Martini <- Hmatrix(A=Amat, G=Gmat, method="Martini", 
#'                      ploidy=4, missingValue=-9, maf=0.05)
#' #Computing H matrix (Munoz)
#' Hmat_Munoz <- Hmatrix(A=Amat, G=Gmat, markers = snp.sol, 
#'                       ploidy=4, method="Munoz",
#'                       roundVar=2,
#'                       missingValue=-9, maf=0.05)
#' }
#'
#' @author Rodrigo R Amadeu, \email{rramadeu@@gmail.com}
#' 
#' @references \emph{Munoz, P. R., Resende, M. F. R., Gezan, S. A., Resende, M. D. V., de los Campos, G., Kirst, M., Huber, D., Peter, G. F. (2014). Unraveling additive from nonadditive effects using genomic relationship matrices. Genetics, 198.4: 1759-1768.}
#' @references \emph{Martini, J. W., Schrauf, M. F., Garcia-Baccino, C. A., Pimentel, E. C., Munilla, S., Rogberg-Munoz, A., ... & Simianer, H. (2018). The effect of the H-1 scaling factors tau and omega on the structure of H in the single-step procedure. Genetics Selection Evolution, 50(1), 16.}
#' 
#' @export

Hmatrix <- function(A=NULL,
                     G=NULL,
                     markers=NULL,
                     c=0,
                     method="Martini",
                     tau=1,
                     omega=1,
                     missingValue=-9,
                     maf=0,
                     ploidy=2,
                     roundVar=3
){
  Aorig <- A
  Gorig <- G
  
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
  
  A <- Aorig#[Anhat,Anhat]
  G <- Gorig[Gnhat,Gnhat]
  
  missingGmatrix <- An[missingGmatrix]
  missingAmatrix <- Gn[missingAmatrix]
  
  Time = as.matrix(proc.time()-Time)
  cat("Completed! Time =", Time[3]/60," minutes \n")
  
  cat("Computing the H matrix... \n")
  
  Time = proc.time()
  
  if(method=="Martini"){
    idA <-rownames(A)
    idG <- rownames(G)
    idH <- unique(c(idG,idA))
    idH <- rev(idH)
    A <- A[idH,idH]
    
    index = is.na(match(idH,idG))
    A11 <- A[index,index]
    A12 <- A[index,!index]
    A21 <- A[!index,index]
    A22 <- A[!index,!index]
    G22 <- G[idH[!index],idH[!index]]
    #if(is.singular.matrix(G22))
    #  stop(deparse("Matrix G22 is singular (not invertible)"))
    A22inv = solve(A22) #A is always invertible
    G22inv = try(solve(G22),silent=TRUE)
    if(inherits(G22inv,"try-error")){
      cat(G22inv)
      stop("G22 not inverting with solve(), try a different/modified G matrix")
    }
    H22 = solve((tau*G22inv+(1-omega)*A22inv))
    H11 = A12 %*% A22inv %*% (H22-A22) %*% A22inv %*% A21  
    H12 = A12 %*% A22inv %*% (H22-A22)
    H21 = (H22-A22) %*% A22inv%*%A21
    H22 = (H22-A22)
    H = A+cbind(rbind(H11,H21),rbind(H12,H22))
    Time = as.matrix(proc.time()-Time)
    cat("\n","Completed! Time =", Time[3]/60," minutes \n")
    return(H)
  }
  
  if(method=="Munoz"){
    A <- Aorig[Anhat,Anhat]
    if(is.null(markers))
      stop("Aborting: For Munoz method you need to specify method object")
    markersmatrix <- Gmatrix(markers,method="MarkersMatrix",ploidy=ploidy,missingValue=missingValue,maf=maf)
    
    #Computing the Variance of G by A classes (A rounded by roundVar)
    classes <- as.numeric(levels(as.factor(A)))
    classes <- unique(round(classes,roundVar))
    n <- length(classes)
    varA <- meanG <- matrix(NA,nrow=nrow(A),ncol=nrow(A))
    varAclasses <- c()
    for(i in 1:n){
      varA[round(A,roundVar)==classes[i]] <- varAclasses[i] <- var(G[round(A,roundVar)==classes[i]])
      meanG[round(A,roundVar)==classes[i]]  <- mean(G[round(A,roundVar)==classes[i]])
    }
    
    varAclasses[varAclasses==0] = NA
    varAclasses[is.infinite(varAclasses)] <- NA
    tmp <- which(is.na(varAclasses))
    for(i in 1:length(tmp)){
      varAclasses[tmp[i]]<-zoo::na.approx(varAclasses)[tmp[i]]
      varA[round(A,roundVar)==classes[tmp[i]]] <- varAclasses[tmp[i]]
    }
    
    #Computaing beta and H
    beta <- 1 - (c+(1/(markersmatrix[Gnhat,Gnhat]))/varA)
    H <- beta*(G-A)+A ######
    Aorig[Anhat,Anhat] = H
    Time = as.matrix(proc.time()-Time)
    cat("\n","Completed! Time =", Time[3]/60," minutes \n")
    cat("\n","Returning H = A matrix corrected by G... \n")
    return(Aorig)
  }
}

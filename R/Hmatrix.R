#########################################################################
# 									#
# Package: AGHmatrix 							#
# 									#
# File: Amatrix.R 							#
# Contains: Amatrix 							#
# 									#
# Written by Rodrigo Rampazo Amadeu 					#
# 									#
# First version: Feb-2014 						#
# Last update: 22-Apr-2015 						#
# License: GNU General Public License version 2 (June, 1991) or later 	#
# 									#
#########################################################################

#' Construction of Relationship Matrix H
#'
#' Given a matrix A and a matrix G returns a H matrix - A matrix corrected by G.
#'
#' @param A A matrix from function Amatrix
#' @param G G matrix from function Gmatrix
#' @param c constant value of H computation
#' @param explore if TRUE performs exploratory analysis of the matrix
#'
#' @return H Matrix with the relationship between the individuals based on pedigree and corrected by molecular information
#'
#' @examples Hmatrix()
#'
#' @author Rodrigo R Amadeu, \email{rramadeu@@gmail.com}
#'
#' @references Patricio Munoz ?
#'
#' @export

Hmatrix <- function(A=NULL,
                    G=NULL,
                    markers=NULL,
                    c=0,
                    explore=FALSE
                    ){
    Aorig <- A
    Gorig <- G

    Time = proc.time()
    cat("Comparing the matrices... \n")
    An <- rownames(Aorig)
    Gn <- rownames(Gorig)
    missingGmatrix <- which(is.na(match(An,Gn)))
    missingAmatrix <- which(is.na(match(Gn,An)))
    Anhat <- An[-missingGmatrix]
    Gnhat <- Gn[-missingAmatrix]
    A <- Aorig[Anhat,Anhat]
    G <- Gorig[Gnhat,Gnhat]

    missingGmatrix <- Gn[missingGmatrix]
    missingAmatrix <- An[missingAmatrix]

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
    beta <- 1 - (c+(1/(markers[Gnhat,Gnhat]))/varA)
    H <- beta*(G-A)+A ######
    Aorig[Anhat,Anhat] = H
    cat("\n","Individuals in A but not in G:",missingAmatrix,"\n")
    cat("\n","Individuals in G but not in A:",missingGmatrix,"\n")
    Time = as.matrix(proc.time()-Time)
    cat("\n","Completed! Time =", Time[3]/60," minutes \n")
    cat("\n","Returning H = A matrix corrected by G...")
    return(Aorig)
}

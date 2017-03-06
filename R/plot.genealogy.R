#########################################################################
# 									
# Package: AGHmatrix 							
# 									
# File: plot.genealogy.R 							
# Contains: plot.genealogy, ancestral, cartesian, cartesian.plot, cartesian.big.plot 
# 									
# Written by Rodrigo Rampazo Amadeu 					
# 									
# First version: 11-Aug-2014 						
# Last update: 06-Mar-2017 						
# License: GNU General Public License version 2 (June, 1991) or later 	
# 									
#########################################################################

#' Construction of Genealogy for an individual
#'
#' Builds a genealogy for a given individual based on his pedigree
#'
#' @param data pedigree data name (data frame in a 3-column way format).
#' @param ind indiviudal name (string)
#' @param gen the max number of generations ago of the plot
#' @param Amatrix relationship matrix in same order of the pedigree (matrix)
#' @param big use it for a better visualization for high gen parameter (TRUE or FALSE)
#' @param fontsize an integer for font size of plot
#' @param angle.names an integer to turn the names of plot  
#' 
#' @return Graphic with ind genealogy until gen generations ago
#'
#' @examples
#' data(ped.mrode)
#' #Build Amatrix diploid (no double reduction proportion)
#' Amat <- Amatrix(data=ped.mrode,ploidy=2,unk=0)
#' #Build Amatrix autotetraploidy (double reduction proportion=0.1)
#' anc.plot(data=ped.mrode,ind="Var4",Amatrix=Amat)
#'
#' @author Rodrigo R Amadeu, \email{rramadeu@@gmail.com}
#'
#' @export

plot.genealogy <- function(data=data,
                           ind=NULL,
                           gen=20,
                           Amat=NULL,
                           fontsize=3,
                           save.image=TRUE,
                           return.plot=TRUE,
                           angle.names=0,
                           big=FALSE){
  
  x.anc <- ancestral(data,ind=ind,gen=gen)
  x<- cartesian.plot(x.anc)
  if(nrow(x) <= 7)
    big <- FALSE
  if(big)
    x<- cartesian.big.plot(x)
  lines.plot <- x
  data.plot <- x[x$data.anc!=0,]
  data.plot$names <- data[,1][data.plot$data.anc]
  if(!is.null(Amatrix)){
    data.plot$inbreed <- round(diag(Amat)[data.plot$data.anc]-1,4)
    inbreed <- TRUE
  }else(inbreed <- FALSE)
  
  p <- ggplot(data.plot, aes(x=xaxis,y=yaxis,label=names)) +
    xlab("Generations Ago") +
    ylab("")
  j<-2
  n <- length(lines.plot$data.anc)
  for(i in 1:n){
    if(!is.na(lines.plot$data.anc[j]))
      if(lines.plot$data.anc[j]!=0){
        p <- p + geom_segment(x=lines.plot$xaxis[i],
                              y=lines.plot$yaxis[i],
                              xend=lines.plot$xaxis[j],
                              yend=lines.plot$yaxis[j],
                              color='gray')
      }
    if(!is.na(lines.plot$data.anc[j+1]))
      if(lines.plot$data.anc[j+1]!=0){
        p <- p +  geom_segment(x=lines.plot$xaxis[i],
                               y=lines.plot$yaxis[i],
                               xend=lines.plot$xaxis[j+1],
                               yend=lines.plot$yaxis[j+1],
                               color='gray')
      }
    j<-j+2
    
  }
  
  p <- p +
    geom_point(colour='gray') +
    geom_text(size=fontsize,angle=angle.names) +
    theme_bw() +
    scale_x_continuous(breaks=1:gen) +
    theme(axis.text.y = element_blank(),
          axis.line.y=element_blank(),
          panel.border=element_blank(),
          panel.grid.major.y=element_blank(),
          panel.grid.minor=element_blank(),
          axis.ticks.y=element_blank())
  if(inbreed)
    p <- p+geom_text(label=data.plot$inbreed,angle=0,vjust=2,size=fontsize,color="#333333")
  
  if(return.plot)
    return(p)
}

ancestral <- function(data,ind,gen=10){
    data <- data.treat(data)
    index <- which(data$ind.data == ind)
    chro <- c(index)
    
    limit=2
    i=1
    while( limit <= gen){
        if(chro[i] == 0){
            chro[length(chro)+1]=0
        }else{
            chro[length(chro)+1] = data$sire[chro[i]]
        }
        
        if(chro[i] == 0){
            chro[length(chro)+1]=0
        }else{
            chro[length(chro)+1] = data$dire[chro[i]]
        }
        
        if( (length(chro)) == -(1-2^limit)){
            if( sum(tail(chro,n=2^(limit-1))) == 0){
                cat("...  \n")
                cat(paste(ind,"has data until generation",limit-1))
                chro <- chro[-(length(chro):(1+length(chro)-2^(limit-1)))]
                return(chro)
            }
            limit=limit+1
        }
        if(limit >= gen){
            print(paste("Reached the limit of generations ",gen))  
            return(chro)            
        }
        i<-i+1
    }
}

#Makes the Cartesian grid
cartesian <- function(gen){
    xaxis <- c(0)
    yaxis <- c(0)
    for(i in 1:gen){
        xaxis <- c(xaxis,rep(i,2^i))
        yaxis <- c(yaxis,rev(seq(2^(i-1):1)),-seq((2^(i-1):1)))
    }
    yaxis <- yaxis
    output<-cbind(xaxis,yaxis)
    return(output)
}

#Costumize the Cartesian grid for larger pedigrees
cartesian.big.plot <- function(x){
    index <- max(x$xaxis)
    x$count <- rep(1,nrow(x))
    big.gen <- big.gen.count <-c()
    for( i in 1:index){
        maxi <- sum(x$data.anc[which(x$xaxis==i)]!=0)
        big.gen.count <- c(big.gen.count,rep(maxi,maxi))
        big.gen <- c(big.gen,seq(1:maxi))
    }
    k=1
    for(i in 2:nrow(x)){
        if(x$data.anc[i] != 0){
            x$yaxis[i] <- big.gen[k]
            x$count[i] <- big.gen.count[k]
            k <- k+1
        }
    }
    x$yaxis <- max(big.gen)*(x$count-x$yaxis)/(x$count-1)
    x$yaxis[1] <- max(big.gen)/2
    x$max <- max(big.gen)[1]
    return(x)
}
    
#Put the data in the cartesian grid        
cartesian.plot <- function(data.anc){
    gen <- log(length(data.anc)+1,base=2)-1
    plane <- cartesian(gen)
    matrix <- cbind(plane,data.anc)
    rownames(matrix) <- c()
    return(data.frame(matrix))
}




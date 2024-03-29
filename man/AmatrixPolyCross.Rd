% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/AmatrixPolyCross.R
\name{AmatrixPolyCross}
\alias{AmatrixPolyCross}
\title{Construction of pedigree-based relationship matrix with parental guessing possibility}
\usage{
AmatrixPolyCross(data = NULL, fixedParent = FALSE)
}
\arguments{
\item{data}{pedigree data name. Unknown value should be equal 0. See example for construction.}

\item{fixedParent}{if false, assumes that all the parents are equally possible parents. If true, assumes that the first parental is known and the others are equally possible parents. Default = FALSE.}
}
\value{
Matrix with the relationship between the individuals.
}
\description{
Creates an additive relationship matrix A based on a non-deterministic pedigree with 4+ columns where each column represents a possible parent. This function was built with the following designs in mind. 
1) A mating design where you have equally possible parents. For example, a generation of insects derived from the mating of three insects in a cage. All the insects in this generation will have the same expected relatedness with all the possible parents (1/3). If there are only two parents in the cage, the function assumes no-inbreeding and the pedigree is deterministic (the individual is offspring of the cross between the two parents). Another example, a population of 10 open-pollinated plants where you harvest the seeds without tracking the mother. 
2) When fixedParent is TRUE: a mating design where you know one parent and might know the other possible parents. For example, a polycross design where you have seeds harvested from a mother plant and possible polen donors.
}
\examples{
#the following pedigree has the id of the individual followed by possible parents
#if 0 is unknown
#the possible parents are filled from left to right
#in the pedigree data frame examples:
#id 1,2,3,4 have unknown parents and are assumed unrelated
#id 5 has three possible parents (1,2,3)
#id 6 has three possible parents (2,3,4)
#id 7 has two parents (deterministic case here, the parents are 3 and 4)
#id 8 has four possible parents (5,6,7,1)

pedigree = data.frame(id=1:8,
                      parent1 = c(0,0,0,0,1,2,3,5),
                      parent2 = c(0,0,0,0,2,3,4,6),
                      parent3 = c(0,0,0,0,3,4,0,7),
                      parent4 = c(0,0,0,0,0,0,0,1),
                      parent5 = 0)

print(pedigree)

AmatrixPolyCross(pedigree)

#when polyCross is set to be true:
#id 5 is offspring of parent 1 in a deterministic way and two other possible parents (2,3)
#id 6 is offspring of parent 2 in a deterministic way and two other possible parents (3,4)
#id 7 has two parents (deterministic case here, the parents are 3 and 4); as before
#id 8 is offspring of parent 5 in a deterministic way and has three other possible parents (6,7,1)

AmatrixPolyCross(pedigree,fixedParent=TRUE)

}
\author{
Rodrigo R Amadeu, \email{rramadeu@gmail.com}
}

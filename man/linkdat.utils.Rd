\name{linkdat.utils}
\alias{spouses}
\alias{offspring}
\alias{descendants}
\alias{pedigreeLoops}
\alias{uninformative}
\title{
linkdat utilities
}
\description{
Utility functions for 'linkdat' objects
}
\usage{
spouses(x, id, original.id)

offspring(x, id, original.id)

descendants(x, id, original.id)

pedigreeLoops(x)

uninformative(x, return.reduced.ped=FALSE)
}
\arguments{
  \item{x}{
a \code{linkdat} object
}
  \item{id}{
an integer, the ID of an individual.
}
 \item{original.id}{
A logical indicating whether 'id' refers to the original ID label or the internal labeling.
}
 \item{return.reduced.ped}{
A logical indicating whether the uninformative individuals should be removed, and the remaining pedigree returned as a new \code{linkdat} object)
}
}
\value{
For \code{spouses(x, id)}, an integer vector containing the ID's of all pedigree members sharing offspring with \code{id}. 

For \code{offspring(x,id)}, a vector containing the ID's of the offspring (i.e. children) of individual \code{id}. 

For \code{descendants(x,id)}, a vector containing the ID's of all descendants (i.e. children, grandchildren, a.s.o.) of individual \code{id}. 

For \code{pedigreeLoops}, a list containing all loops found in the pedigree. Each loop is represented as a list with elements 'top', a 'bottom' individual, 'pathA' (individuals forming a path from top to bottom) and 'pathB' (creating a different path from top to bottom, with no individuals in common with pathA). 

For \code{uninformative}, either a numeric vector (if \code{return.reduced.ped=FALSE}) containing the ID labels of the uninformative individuals, or (if \code{return.reduced.ped=TRUE}) a \code{linkdat} object where the uninformative individuals have been removed.
} 
\author{
Magnus Dehli Vigeland
}

\examples{
p=cbind(ID=2:9, FID=c(0,0,2,0,4,4,0,2), MID=c(0,0,3,0,5,5,0,8), 
        SEX=c(1,2,1,2,1,2,2,2), AFF=c(2,1,2,1,2,1,1,2))
x=linkdat(p)

desc = descendants(x, 2)
spous = spouses(x, 2)
offs = offspring(x, 2)
stopifnot(setequal(spous, c(3,8)), setequal(offs, c(4,9)), setequal(desc, c(4,6,7,9)))

#Creating a loop and detecting it with 'pedigreeLoops' 
#(note that we get two loops, one for each inbred child):
loopx = addOffspring(x, father=4, mother=9, noffs=2)
lps = pedigreeLoops(loopx)
stopifnot(lps[[1]]$top==2, setequal(sapply(lps, '[[', 'bottom'), 10:11))

#an example with uninformative individuals
data(dominant)
x = linkdat(dominant)
un = uninformative(x)
stopifnot(setequal(un, c(21,23)))
}
\keyword{ math }
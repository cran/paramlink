\name{linkdat.utils}
\alias{spouses}
\alias{offspring}
\alias{descendents}
\alias{ped.dist.matrix}
\alias{getpaths}
\title{
linkdat utilities
}
\description{
Utility functions for linkdat objects
}
\usage{
spouses(x, i)

offspring(x, i)

descendents(x, i)

ped.dist.matrix(x)

getpaths(obj, i,j)
}
\arguments{
  \item{x}{
a \code{linkdat} object
}
  \item{i, j}{
an integer, the ID of an individual.
}
  \item{obj}{
a matrix of pairwise path distances, or a \code{linkdat} object.
}
}
\value{
For \code{spouses(x,i)}, an integer vector contining the ID's of all pedigree members sharing offspring with \code{i}. 

For \code{offspring(x,i)}, a vector containing the ID's of the offspring (i.e. children) of individual i. 

For \code{descendents(x,i)}, a vector containing the ID's of all descendents (i.e. children, grandchildren, a.s.o.) of individual i. 

For \code{ped.dist.matrix(x)}, a symmetrix integer matrix (d_ij), where d_ij is the pedigree distance (i.e., the length of the shortest pedigree path) between i and j. For example, pedigree members with the relations below would have the indicated pedigree distances:
parent/offspring	1,
siblings			2
parents				2
uncle/niece			3
first cousins		4

For \code{getpaths(x,i,j)}, a data.frame with \code{ped.dist.matrix(x)[i,j]+1} columns. Each row is a shortest pedigree path (endpoints included) from \code{i} to \code{j}. 
} 
\author{
Magnus Dehli Vigeland
}

\seealso{
\code{\link{subset.linkdat}}
}
\examples{
p=cbind(ID=1:6, FID=c(0,0,1,0,3,3), MID=c(0,0,2,0,4,4), 
        SEX=c(1,2,1,2,1,2), AFF=c(2,1,2,1,2,1))
x=linkdat(p, model=1)

#all descendents of the grandfather:
descendents(x,1)

#extracting the nuclear family where individual 3 is the father:
nucl=c(3, spouses(x,3), offspring(x,3))
y=subset(x, subset=nucl)
y
summary(y)

d=ped.dist.matrix(x)
d
#There should be 2 paths of length 3 from individual 1 to individual 4:
paths=getpaths(d,1,4)
stopifnot(dim(paths)==c(2,4))

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ math }

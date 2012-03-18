\name{markers}
\alias{setMarkers}
\alias{getMarkers}
\alias{modifyMarker}
\alias{mendelianCheck}

\title{
Marker genotypes for linkdat objects  
}
\description{
Functions for setting and manipulating marker genotypes for 'linkdat' objects.
}
\usage{
setMarkers(x, m, missing = 0)
getMarkers(x, markers)
modifyMarker(x, id, alleles, marker=1)
mendelianCheck(x)
}
\arguments{
  \item{x}{
a \code{linkdat} object
}
  \item{m}{
a matrix with one row per pedigree member and two columns per marker, containing the marker alleles. 
}
  \item{missing}{
a numeric - or character - of length 1, indicating the code for missing alleles.
}
  \item{markers}{
an integer vector containing marker IDs.
}
  \item{id}{
an integer vector, indicating the individual(s) to be modified.
}
  \item{alleles}{
a vector of length 2 (if \code{id} is a single integer), or a matrix with \code{length(id)} rows and 2 columns.
}
  \item{marker}{
an integer indicating which marker to modify.
}
}
\value{
For \code{setMarkers} and \code{modifyMarkers}, a \code{linkdat} object is returned, whose \code{markerdata} element has been set/modified. 
The \code{markerdata} contains the marker allele matrix (numerically coded) and has the following two attributes: 
\item{alleles}{The allele names. At the moment the notation for these must be common for all markers.}
\item{missing}{The code for missing data.}

For \code{getMarkers}, a matrix with \code{x$nInd} rows and \code{2*length(markers)} columns.
For \code{mendelianCheck}, a numeric containing the indices of the markers that did not pass the tests.
}
\author{
Magnus Dehli Vigeland}

\seealso{
\code{\link{linkdat}}
}
\examples{
data(Xped)
x=linkdat(Xped, model=4)
#setting 2 markers, both homozygous for all individuals 
x=setMarkers(x, m=matrix(rep(1:2, each=30), nrow=15))
x

#introducing a some mendelian errors:
x = modifyMarker(x, id=13, alleles=c(1,2), marker=1)
x = modifyMarker(x, id=14, alleles=c(2,2), marker=1)
x = modifyMarker(x, id=15, alleles=c(1,2), marker=1)
plot(x, marker=1)
mendelianCheck(x)
}
\keyword{ manip }

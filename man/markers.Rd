\name{markers}
\alias{marker}
\alias{addMarker}
\alias{setMarkers}
\alias{getMarkers}
\alias{removeMarkers}
\alias{modifyMarker}
\alias{modifyMarkerMatrix}
\alias{swapGenotypes}
\alias{mendelianCheck}
\alias{setMap}

\title{
Marker functions  
}
\description{
Functions for setting and manipulating marker genotypes for 'linkdat' objects.
}
\usage{
marker(x, ..., allelematrix, alleles=NULL, afreq=NULL, missing=0,
       chrom=NA, name=NA, pos=NA)

addMarker(x, m, alleles=NULL, afreq=NULL, missing=0)

setMarkers(x, m, annotations=NULL, missing=0)

getMarkers(x, markernames=NULL, chroms=NULL, from=NULL, to=NULL)

removeMarkers(x, markers, markernames=NULL, chroms=NULL, from=NULL, to=NULL)

modifyMarker(x, marker, ids, genotype, alleles, afreq, chrom, name, pos)

modifyMarkerMatrix(x, ids, new.alleles)

swapGenotypes(x, ids)

mendelianCheck(x, remove=FALSE, verbose=!remove)
}
\arguments{
  \item{x}{
a \code{\link{linkdat}} object
}
  \item{...}{
an even number of vectors, indicating individuals and their genotypes. See examples.
}
  \item{allelematrix}{
a matrix with one row per pedigree member and two columns per marker, containing the alleles for a single marker. 
}
  \item{m}{
a \code{marker} object or a matrix with alleles. (In \code{setMarkers} this matrix can contain data of several markers.)
}
  \item{missing}{
a numeric - or character - of length 1, indicating the code for missing alleles.
}
  \item{chrom}{
NA or an integer (the chromosome number of the marker).
}
  \item{name}{
NA or a character (the name of the marker).
}
  \item{pos}{
NA or a non-negative real number indicating the genetic position (in cM) of the marker.
}
  \item{annotations}{
a list of marker annotations.
}
  \item{markernames}{
NULL or a character vector.
}
  \item{chroms}{
NULL or a numeric vector of chromosome numbers.
}
  \item{from}{
NULL or a single numeric.
}
  \item{to}{
NULL or a single numeric.
}
  \item{marker, markers}{
a numeric indicating which marker(s) to use/modify.
}
  \item{ids}{
a numeric indicating individual(s) to be modified. In \code{swapGenotypes} this must have length 2.
}
  \item{genotype}{
a vector of length 1 or 2, containing the genotype to be given the \code{ids} individuals. See examples.
}
  \item{alleles}{
a numeric or character vector containing allele names.
}
  \item{afreq}{
a numerical vector with allele frequencies. An error is given if they don't sum to 1 (rounded to 3 decimals).
}
 \item{new.alleles}{
a numerical matrix of dimensions \code{length(ids), 2*x$nMark}. Entries refer to the internal allele numbering.
}
 \item{remove}{
a logical. If FALSE, the function returns the indices of markers found to incorrect. 
If TRUE, a new \code{linkdat} object is returned, where the incorrect markers have been deleted.
}
 \item{verbose}{
a logical. If TRUE, details of the markers failing the tests are shown.
}
}
\value{
The \code{marker} function returns an object of class \code{marker}: This is a numerical 2-column matrix with one row per individual, and attributes 
'alleles' (a character vector with allele names), 'nalleles' (the number of alleles) and 'missing' (the input symbol for missing marker alleles), 'chrom' (chromosome number), 'name' (marker identifier), 'pos' (chromosome position in cM).

For \code{addMarker}, \code{setMarkers}, \code{removeMarkers}, \code{modifyMarker}, \code{modifyMarkerMatrix} and \code{swapGenotypes}, a \code{linkdat} object is returned, whose \code{markerdata} element has been set/modified. 

For \code{getMarkers} a numeric vector containing marker numbers (i.e. their indices in \code{x$markerdata}) for the markers whose
"name" attribute is contained in \code{markernames}, "chrom" attribute is contained in \code{chroms}, and "pos" attribute is between
\code{from} and \code{to}. NULL arguments are skipped, so \code{getMarkers(x)} will return \code{seq_len(x$nMark)} (i.e. all markers).

For \code{mendelianCheck}, a numeric containing the indices of the markers that did not pass the tests, or (if \code{remove=TRUE}) a new \code{linkdat} object where the failing markers are removed.
}
\author{
Magnus Dehli Vigeland}

\seealso{
\code{\link{linkdat}}
}
\examples{
x = linkdat(toyped)
x = removeMarkers(x,1) # removing the only marker of x.
x

# creating a SNP with alleles 'a' and 'b', for which individual 1 is heterozygous, 
# 2 and 4 are homozygous for the 'b' allele, and missing genotype for individual 3:
m1 = marker(x, 1, c('a','b'), c(2,4), c('b','b')) 
x = addMarker(x, m1)

# a rare SNP for which all individuals are homozygous for the common allele.
# The 'alleles' argument must be specified since only one allele is present.
m2 = marker(x, 1:4, 1, alleles=1:2, afreq=c(0.99, 0.01))
x = addMarker(x, m2)
x

# creating a Mendelian error!
x = modifyMarker(x, marker=1, ids=3, genotype='a') 

# making the second marker triallelic, with heterozygous parents.
x = modifyMarker(x, marker=2, ids=1:2, genotype=1:2, alleles=1:3)
x

plot(x, marker=1)
err_index = mendelianCheck(x, remove=FALSE)
xnew = mendelianCheck(x, remove=TRUE)
stopifnot(err_index==1, all.equal(xnew, removeMarkers(x, err_index)))
}
\keyword{ manip }

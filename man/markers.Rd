\name{markers}
\alias{marker}
\alias{addMarker}
\alias{setMarkers}
\alias{removeMarkers}
\alias{modifyMarker}
\alias{mendelianCheck}
\alias{setMap}

\title{
Marker genotypes for linkdat objects  
}
\description{
Functions for setting and manipulating marker genotypes for 'linkdat' objects.
}
\usage{
marker(x, ..., allelematrix, alleles=NULL, afreq=NULL, missing=0,
       chrom=NA, name=NA, pos=NA)
addMarker(x, m, alleles=NULL, afreq=NULL, missing=0)
setMarkers(x, m, map=NULL, dat=NULL, freq=NULL, missing = 0)
removeMarkers(x, markers)
modifyMarker(x, marker, ids, genotype, alleles, afreq, chrom, name, pos)
mendelianCheck(x, verbose=TRUE)
}
\arguments{
  \item{x}{
a \code{linkdat} object
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
  \item{marker, markers}{
a numeric indicating which marker(s) to use/modify.
}
  \item{ids}{
a numeric indicating individual(s) to be modified.
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
  \item{map}{
a character containing the path to a map file in MERLIN style, i.e. with three columns: CHR (chromosome number), MARKER (marker identifier) and POS (centiMorgan position of the marker).
}
 \item{dat}{
a character containing the path to a dat file, MERLIN style. See MERLIN's online tutorial for details on format. Our use of this file here is mainly to communicate the ordering of the markers.
}
 \item{freq}{
a character with the path to a allele frequency file in MERLIN (short) format. (Not used if \code{map=NULL}.)
}
 \item{verbose}{
a logical. If TRUE, details of the markers failing the tests are shown.
}
}
\value{
The \code{marker} function returns an object of class \code{marker}: This is a numerical 2-column matrix with one row per individual, and attributes 
'alleles' (a character vector with allele names), 'nalleles' (the number of alleles) and 'missing' (the input symbol for missing marker alleles), 'chrom' (chromosome number), 'name' (marker identidier), 'pos' (chromosome position in cM).

For \code{addMarker}, \code{setMarkers}, \code{removeMarkers} and \code{modifyMarkers}, a \code{linkdat} object is returned, whose \code{markerdata} element has been set/modified. 

For \code{mendelianCheck}, a numeric containing the indices of the markers that did not pass the tests.
}
\author{
Magnus Dehli Vigeland}

\seealso{
\code{\link{linkdat}}
}
\examples{
data(toyped)
x = linkdat(toyped)
x
x = removeMarkers(x,1) # removing the only marker of x.

# creating a SNP with alleles 'a' and 'b', for which individual 1 is heterozygous, 2 and 4 are homozygous for the 'b' allele, 
# and missing genotype for individual 3:
m1 = marker(x, 1, c('a','b'), c(2,4), c('b','b')) 
m1
x = addMarker(x, m1)

# a rare SNP for which all individuals are homozygous for the common allele.
m2 = marker(x, 1:4, 1, alleles=1:2, afreq=c(0.99, 0.01)) # the 'alleles' argument is important since only one allele is present.
m2
x = addMarker(x, m2)
x

x = modifyMarker(x, marker=1, ids=3, genotype='a') # this is a Mendelian error!
x = modifyMarker(x, marker=2, ids=1:2, genotype=1:2, alleles=1:3) # making the second marker triallelic, with heterozygous parents.
x

plot(x, marker=1)
err = mendelianCheck(x)
stopifnot(err==1)
x = removeMarkers(x, err)
}
\keyword{ manip }

\name{linkdat}
\alias{linkdat}
\alias{print.linkdat}
\alias{summary.linkdat}
\alias{as.data.frame.linkdat}
\alias{subset.linkdat}
\alias{write.linkdat}
\title{
linkdat objects
}
\description{
Functions to create and display 'linkdat' objects.
}
\usage{
linkdat(ped, model=NULL, map=NULL, dat=NULL, freq=NULL, verbose=TRUE, 
        missing=0, ...)

\method{print}{linkdat}(x, ..., markers)

\method{summary}{linkdat}(object, ...)

\method{as.data.frame}{linkdat}(x, ..., famid=F, markers=seq_len(x$nMark), 
        alleles=NULL, missing=NULL, singleCol=FALSE, sep="")

\method{subset}{linkdat}(x, subset=x$orig.ids, ..., markers=seq_len(x$nMark))

write.linkdat(x, prefix="", what=c("ped", "map", "dat", "freq", "model"), 
              merlin=FALSE)
}
\arguments{
 \item{ped}{
a matrix, data.frame or a character with the path to a pedigree file in standard LINKAGE format.
}
 \item{model}{
either a \code{linkdat.model} object (typically \code{y$model} for some linkdat object \code{y}), or a single integer with the following meaning: 
1 = autosomal dominant; 2 = autosomal recessive; 3 = X-linked dominant; 4 = X-linked recessive.
In each of these cases, the disease is assumed fully penetrant and the disease allele frequency is set to 0.00001. If \code{model=NULL}, no model is set.
}
 \item{map}{
a character with the path to a map file in MERLIN format. If non-NULL, a dat file must also be given (next item).
}
 \item{dat}{
a character with the path to a dat file in MERLIN format. (Only needed if \code{map} is non-NULL.)
}
 \item{freq}{
a character with the path to a allele frequency file in MERLIN (short) format. If NULL, all markers are interpreted as equifrequent.
}
 \item{verbose}{
a logical: verbose output or not.
}
 \item{markers}{
a numeric indicating which markers should be included/printed.
}
  \item{x, object}{
a \code{linkdat} object.
}
  \item{famid}{
a logical indicating if the family identifier should be included as the first column.
}
  \item{alleles}{
a character containing the allele names, e.g. \code{alleles=c("A","B")}.
}
  \item{missing}{
a character of length 1 indicating how missing alleles should be represented.
}
  \item{singleCol}{
a logical: Should the two alleles for each marker be pasted into one column or kept in separate columns?
}
  \item{sep}{
a single character to be used as allele separator if \code{singleCol=TRUE}.
}
  \item{subset}{
a numeric containing the individuals in the sub-pedigree to be extracted. NB: No pedigree checking is done here, so make sure the subset form a meaningful, closed pedigree.
}
  \item{prefix}{
a character string giving the prefix of the files. For instance, if \code{prefix="fam1"} and \code{what=c("ped", "map")}, the files "fam1.ped" and "fam1.map" will be created.
}
  \item{what}{
a character vector forming a subset of c("ped", "map", "dat", "freq", "model"), indicating which files should be created. All files are written in MERLIN style (but see the next item!)
}
  \item{merlin}{
a logical. If TRUE, the marker alleles are relabeled to 1,2,..., making sure that the generated files are readable by MERLIN (which does not accept non-numerical allele labels in the frequency file.) If FALSE (the default)
the allele labels are unchanged. In this case, \code{x} should be exactly reproducible from the files. (See examples.)
}
  \item{\dots}{
further arguments.
}
} 
\details{
The file (or matrix or data.frame) \code{ped} must describe one or several pedigrees in standard LINKAGE format, i.e. with the following columns in correct order:

 1 Family id (optional) (FAMID)
 
 2 Individual id (ID),
 
 3 Father id (FID),
 
 4 Mother id (MID),
 
 5 Gender (SEX): 1 = male, 2 = female,

 6 Affection status (AFF): 1 = unaffected, 2 = affected, 0 = unknown,

 7 First allele of first marker,

 8 Second allele of first marker,

 9 First allele of second marker,

 a.s.o.

Only columns 2-6 are mandatory. The first column is automatically interpreted as family id if it has repeated elements.

Internally the individuals are relabeled as 1,2,..., but this should rarely be of concern to the end user. Some pedigree checking is done, but it is recommended to plot the pedigree before doing any analysis.

Details on the formats of map, dat and frequency files can be found in the online MERLIN tutorial: \url{http://www.sph.umich.edu/csg/abecasis/Merlin/}
}

\value{
A \code{linkdat} object, or a list of \code{linkdat} objects. A linkdat object is essentially a list with the following entries, some of which can be NULL.
 \item{pedigree }{\code{data.frame} with 5 columns (ID, FID, MID, SEX, AFF) describing the pedigree in linkage format. (NB: Internal labeling used.)}
 \item{orig.ids}{the original individual id labels.}
 \item{nInd}{the number of individuals in the pedigree.}
 \item{founders}{vector of the founder individuals. (NB: Internal labeling used.)}
 \item{nonfounders}{vector of the nonfounder individuals (NB: Internal labeling used.)}
 \item{subnucs}{list containing all (maximal) nuclear families in the pedigree. Each nuclear family is given as a vector of the form c(pivot, father, mother, child1, ...), where the pivot is either the id of the individual linking the nuclear family to the rest of the pedigree, or 0 if there are none. (NB: Internal labeling used.)}
 \item{markerdata}{a list of \code{\link{marker}} objects.}
 \item{nMark}{the number of markers.}
 \item{available}{a numeric vector containing IDs of available individuals. Used for simulations and plots.}
 \item{model}{a \code{linkdat.model} object, essentially a list containing the model parameters. See \code{\link{setModel}} for details.}
 \item{initial_probs}{(Not relevant for end users.)}
}

\author{
Magnus Dehli Vigeland
}

\seealso{
\code{\link{setModel}}, \code{\link{modifyPedigree}}, \code{\link{linkdat.utils}}
}
\examples{
data(toyped)
toyped
x = linkdat(toyped, model=1)
x
summary(x)

# test read/write:
x = modifyMarker(x, 1, alleles=c('B','C'), afreq=c(.9, .1), chrom=2, name="SNP1", pos=123)
write.linkdat(x, prefix="toy")
y = linkdat("toy.ped", map="toy.map", dat="toy.dat", freq="toy.freq", model=1)
unlink(c("toy.ped", "toy.map", "toy.dat", "toy.freq", "toy.model"))
stopifnot(isTRUE(all.equal(x,y)))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ IO }
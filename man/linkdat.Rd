\name{linkdat}
\alias{linkdat}
\alias{singleton}
\alias{print.linkdat}
\alias{summary.linkdat}
\alias{subset.linkdat}
\alias{write.linkdat}
\title{
linkdat objects
}
\description{
Functions to create and display 'linkdat' objects.
}
\usage{
linkdat(ped, model=NULL, map=NULL, dat=NULL, freq=NULL, annotations=NULL, 
        missing=0, header=FALSE, checkped=TRUE, verbose=TRUE, ...)
        
singleton(id, sex=1, famid=1, verbose=FALSE, ...)

\method{print}{linkdat}(x, ..., markers)

\method{summary}{linkdat}(object, ...)

\method{subset}{linkdat}(x, subset=x$orig.ids, ..., markers=seq_len(x$nMark))

write.linkdat(x, prefix="", what=c("ped", "map", "dat", "freq", "model"), 
              merlin=FALSE)
}
\arguments{
 \item{ped}{
a matrix, data.frame or a character with the path to a pedigree file in standard LINKAGE format. (See details)
}
 \item{model}{
either a \code{linkdat.model} object (typically \code{y$model} for some linkdat object \code{y}), or a single integer with the following meaning: 
1 = autosomal dominant; 2 = autosomal recessive; 3 = X-linked dominant; 4 = X-linked recessive.
In each of these cases, the disease is assumed fully penetrant and the disease allele frequency is set to 0.00001. If \code{model=NULL}, no model is set.
}
 \item{map}{
a character with the path to a map file in MERLIN format, or NULL. If non-NULL, a dat file must also be given (next item).
}
 \item{dat}{
a character with the path to a dat file in MERLIN format, or NULL. (Only needed if \code{map} is non-NULL.)
}
 \item{freq}{
a character with the path to a allele frequency file in MERLIN (short) format, or NULL. If NULL, all markers are interpreted as equifrequent.
}
 \item{annotations}{
a list (of the same length and in the same order as the marker columns in \code{x}) of marker annotations. If this is non-NULL, then all of \code{map, dat, freq} should be NULL.
}
 \item{missing}{
the character (of length 1) used for missing alleles. Defaults to "0".
}
 \item{header}{
a logical, relevant only if \code{ped} points to a ped file: If TRUE, the first line of the ped file is skipped.
}
 \item{checkped}{
a logical. If FALSE, no checks for pedigree errors are performed.
}
 \item{verbose}{
a logical: verbose output or not.
}
 \item{id, sex}{
single numerics describing the individual ID and gender of the singleton.
}
 \item{markers}{
a numeric indicating which markers should be included/printed.
}
  \item{x, object}{
a \code{linkdat} object.
}
  \item{famid}{
a numeric: the family ID of the singleton. 
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

Details on the formats of map, dat and frequency files can be found in the online MERLIN tutorial: \url{http://csg.sph.umich.edu/abecasis/Merlin/}

A singleton is a special \code{linkdat} object whose pedigree contains 1 individual. The class attribute of a singleton is \code{c('singleton', 'linkdat')}
}

\value{
A \code{linkdat} object, or a list of \code{linkdat} objects. A linkdat object is essentially a list with the following entries, some of which can be NULL.
 \item{pedigree }{\code{data.frame} with 5 columns (ID, FID, MID, SEX, AFF) describing the pedigree in linkage format. (NB: Internal labeling used.)}
 \item{orig.ids}{the original individual id labels.}
 \item{nInd}{the number of individuals in the pedigree.}
 \item{founders}{vector of the founder individuals. (NB: Internal labeling used.)}
 \item{nonfounders}{vector of the nonfounder individuals (NB: Internal labeling used.)}
 \item{hasLoops}{a logical: TRUE if the pedigree is inbred.}
 \item{subnucs}{list containing all (maximal) nuclear families in the pedigree. Each nuclear family is given as a vector of the form c(pivot, father, mother, child1, ...), where the pivot is either the id of the individual linking the nuclear family to the rest of the pedigree, or 0 if there are none. (NB: Internal labeling used.)}
 \item{markerdata}{a list of \code{\link{marker}} objects.}
 \item{nMark}{the number of markers.}
 \item{available}{a numeric vector containing IDs of available individuals. Used for simulations and plots.}
 \item{model}{a \code{linkdat.model} object, essentially a list containing the model parameters. See \code{\link{setModel}} for details.}
 \item{loop_breakers}{a matrix with original loop breaker ID's in the first column and their duplicates in the second column. This is set by \code{\link{breakLoops}}.}
}

\author{
Magnus Dehli Vigeland
}

\seealso{
\code{\link{setModel}}, \code{\link{modifyPedigree}}, \code{\link{linkdat.utils}}
}
\examples{
x = linkdat(toyped, model=1)
x
summary(x)

#### test read/write:
x = modifyMarker(x, 1, alleles=c('B','C'), afreq=c(.9, .1), chrom=2, name="SNP1", pos=123)
write.linkdat(x, prefix="toy")
y = linkdat("toy.ped", map="toy.map", dat="toy.dat", freq="toy.freq", model=1)
unlink(c("toy.ped", "toy.map", "toy.dat", "toy.freq", "toy.model"))
stopifnot(isTRUE(all.equal(x,y)))

#### test singletons:
w = singleton(id=3, sex=2)
T1 = all.equal(w, linkdat(ped=rbind(c(3,0,0,2,1))))
w = markerSim(w, N=5, alleles=2, afreq=c(0.1,.9))
T2 = all.equal(w, relabel(relabel(w, 10), 3))
T3 = all.equal(w, swapSex(swapSex(w, 3), 3))
T4 = all.equal(w, swapAff(swapAff(w, 3), 3))
stopifnot(T1, T2, T3, T4)

#### several ways of creating the same linkdat object:
alleles = c(157,160,163)
afreq = c(0.3, 0.3, 0.4)
gt10 = c(160, 160)
gt14 = c(160, 163)

z1 = relabel(addOffspring(nuclearPed(1), father=3, noffs=1, aff=2), 10:14)
z1 = addMarker(z1, marker(z1, 10, gt10, 14, gt14, alleles=alleles, afreq=afreq))
z1 = setModel(z1, 2)

z2 = addParents(relabel(nuclearPed(1), 12:14), 12, father=10, mother=11)
z2 = addMarker(z2, rbind(gt10, 0, 0, 0, gt14), alleles=alleles, afreq=afreq)
z2 = setModel(swapAff(z2, 14), 2)

z3 = linkdat(data.frame(ID=10:14, FID=c(0,0,10,0,12), MID=c(0,0,11,0,13), 
             SEX=c(1,2,1,2,1), AFF=c(1,1,1,1,2), 
             M=c("160/160", "0/0", "0/0", "0/0", "160/163")), model=2)
z3 = modifyMarker(z3, 1, alleles=alleles, afreq=afreq)

write.linkdat(z1, prefix="test")
z4 = linkdat("test.ped", map="test.map", dat="test.dat", freq="test.freq", 
             model=2)
z4 = modifyMarker(z4, 1, alleles=alleles, chrom=NA, pos=NA, name=NA)

write.linkdat(z1, prefix="test", merlin=TRUE)
z5 = linkdat("test.ped", map="test.map", dat="test.dat", freq="test.freq", 
             model=2)
z5 = modifyMarker(z5, 1, alleles=alleles, chrom=NA, pos=NA, name=NA)

stopifnot(isTRUE(all.equal(z1,z2)), isTRUE(all.equal(z1,z3)), 
          isTRUE(all.equal(z1,z4)), isTRUE(all.equal(z1,z5)))
unlink(c("test.ped", "test.map", "test.dat", "test.freq", "test.model"))
}

\keyword{ IO }
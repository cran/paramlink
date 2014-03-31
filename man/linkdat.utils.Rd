\name{linkdat.utils}
\alias{linkdat.utils}
\alias{ancestors}
\alias{descendants}
\alias{cousins}
\alias{grandparents}
\alias{nephews_nieces}
\alias{offspring}
\alias{parents}
\alias{siblings}
\alias{spouses}
\alias{unrelated}
\alias{related.pairs}
\alias{pedigreeLoops}
\alias{breakLoops}
\alias{tieLoops}
\alias{all.equal.linkdat}
\title{
linkdat utilities
}
\description{
Utility functions for 'linkdat' objects
}
\usage{
ancestors(x, id)
descendants(x, id, original.id=TRUE)

cousins(x, id, degree=1, removal=0, half=NA, original.id=TRUE)
grandparents(x, id, degree=2, original.id=TRUE)
nephews_nieces(x, id, removal=1, half=NA, original.id=TRUE)
parents(x, id, original.id=TRUE)
offspring(x, id, original.id=TRUE)
siblings(x, id, half=NA, original.id=TRUE)
spouses(x, id, original.id=TRUE)
unrelated(x, id, original.id=TRUE)

related.pairs(x, relation=c('parents', 'siblings', 'grandparents', 'nephews_nieces', 
              'cousins', 'spouses', 'unrelated'), available=F, original.id=T, 
              interfam=c("none", "founders", "all"), ...)

pedigreeLoops(x)
breakLoops(x, loop_breakers)
tieLoops(x)

\method{all.equal}{linkdat}(target, current, ...)
}
\arguments{
  \item{x}{
a \code{\link{linkdat}} object. In \code{related.pairs} possibly a list of \code{linkdat} objects.
}
  \item{id}{
a numerical ID label.
}
 \item{original.id}{
a logical indicating whether 'id' refers to the original ID label or the internal labeling.
}
 \item{degree}{
a non-negative integer.
}
 \item{removal}{
a non-negative integer
}
 \item{half}{
a logical or NA. If TRUE (resp FALSE), only half (resp. full) siblings/cousins/nephews/nieces are returned. If NA, both categories are inclucded. 
}
 \item{relation}{
one of the words (possibly truncated) \code{parents}, \code{siblings}, \code{grandparents}, \code{nephews_nieces}, \code{cousins}, \code{spouses}, \code{unrelated}.
}
\item{interfam}{
one of the words (possibly truncated) \code{none}, \code{founders} or \code{all}, specifying which interfamiliar pairs should be included as unrelated in 
the case where \code{x} is a list of several pedigrees. If \code{none}, only intrafamiliar pairs are considered; 
if \code{founders} all interfamiliar pairs of (available) founders are included; if \code{all}, all interfamiliar (available) pairs are included.
}
 \item{available}{
a logical, if TRUE only pairs of available individuals are returned. 
}

 \item{loop_breakers}{
A vector with IDs of the loop breaker individuals.
}
 \item{target, current}{
Two \code{linkdat} objects to be compared.
}
 \item{...}{
Not used.
}
}
\value{
For \code{ancestors(x,id)}, a vector containing the ID's of all ancestors of the individual \code{id}. 
For \code{descendants(x,id)}, a vector containing the ID's of all descendants (i.e. children, grandchildren, a.s.o.) of individual \code{id}. 

The functions \code{cousins}, \code{grandparents}, \code{nephews_nieces}, \code{offspring}, \code{parents}, \code{siblings}, \code{spouses}, \code{unrelated}, 
each returns an integer vector containing the ID's of all pedigree members having the specified relationship with \code{id}. 

For \code{related.pairs} a matrix with two columns. Each row gives of a pair of pedigree members with the specified relation. If the input is a list of multiple
pedigrees, the matrix entries are characters of the form "X-Y" where X is the family ID and Y the individual ID of the person.

For \code{pedigreeLoops}, a list containing all loops found in the pedigree. 
Each loop is represented as a list with elements 'top', a 'bottom' individual, 'pathA' (individuals forming a path from top to bottom) and 'pathB' 
(creating a different path from top to bottom, with no individuals in common with pathA). Note that the number of loops reported here
counts all closed paths in the pedigree and will in general be larger than the genus of the underlying graph.

For \code{breakLoops}, a \code{linkdat} object in which the indicated loop breakers are duplicated. The returned object will also have a non-null \code{loop_breakers} entry, namely a matrix with the IDs of the original loop breakers in the first column and the duplicates in the second.

For \code{tieLoops}, a \code{linkdat} object in which any duplicated individuals (as given in the \code{x$loop_breakers} entry) are merged. Ideally, \code{tieLoops(breakLoops(x, loop_breakers))} should return \code{x} for any linkdat object \code{x}.
} 
\author{
Magnus Dehli Vigeland
}

\examples{
p = cbind(ID=2:9, FID=c(0,0,2,0,4,4,0,2), MID=c(0,0,3,0,5,5,0,8), 
        SEX=c(1,2,1,2,1,2,2,2), AFF=c(2,1,2,1,2,1,1,2))
x = linkdat(p)

desc = descendants(x, 2)
spous = spouses(x, 2)
offs = offspring(x, 2)

stopifnot(setequal(spous, c(3,8)), setequal(offs, c(4,9)), setequal(desc, c(4,6,7,9)))

# Creating a loop and detecting it with 'pedigreeLoops' 
# (note that we get two loops, one for each inbred child):
loopx = addOffspring(x, father=4, mother=9, noffs=2)
lps = pedigreeLoops(loopx)
stopifnot(lps[[1]]$top == 2, setequal(sapply(lps, '[[', 'bottom'), 10:11))

# We add genotypes for a single SNP marker and compute a LOD score under a dominant model.
loopx = setMarkers(loopx, cbind(1,c(2,1,2,1,2,1,1,2,1,1)))
loopx = setModel(loopx, 1)
LOD1 = lod(loopx, theta=0.1, loop_breaker=9)
stopifnot(round(LOD1, 3) == 1.746)

# We could also break the loop manually before computing the LOD:
loopfree = breakLoops(loopx, 9)
LOD2 = lod(loopfree, theta=0.1)
stopifnot(all.equal(loopx, tieLoops(loopfree)))
stopifnot(all.equal(LOD1, LOD2))
}
\keyword{ math }
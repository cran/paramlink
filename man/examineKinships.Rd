\name{examineKinships}
\alias{examineKinships}
\alias{IBDestimate}

\title{
Relatedness estimation
}
\description{
The function \code{IBDestimate} estimates IBD coefficients (k0, k1, k2) for specified pairs of pedigree members, 
using maximum likelihood methods (imported from the maxLik package). 

Checking for pedigree errors is conveniently done with the \code{examineKinships} function: This computes IBD 
coefficients for all pairs in the pedigree(s) and produces a color-coded plot where wrong relationships 
are easy to spot.

These functions work with SNP markers only.
}
\usage{
examineKinships(x, who = "all", interfam = c("founders", "none", "all"), 
                makeplot = T, ...)

IBDestimate(g1, g2=NULL, dat, f=NULL, error=NULL, plot.action=2, pointcol=2, cex=1, ...)
}
\arguments{
  \item{x}{
A \code{\link{linkdat}} object, or a list of such.
}
  \item{who}{
A character vector of one or more of the words "parents", "siblings", "grandparents", "cousins", 
"distant" and "unrelated". Two short forms are possible: "all" (all of the above) and "close" 
(all of the above except "distant" and "unrelated".)
}
  \item{interfam}{
A character; either "founders", "none" or "all", indicating which interfamiliar pairs of individuals
should be included. Only relevant if \code{x} is a list of several \code{linkdat} objects.
}
  \item{makeplot}{
A logical.
}
  \item{g1, g2}{
The pairs of individuals for which IBD estimates should be computed. Several formats are possible. 
If \code{g1} is a vector of ID labels, and \code{g2} NULL, then all pairs among the \code{g1} are 
included. If \code{g2} is not NULL, only pairs intersecting both \code{g1} and \code{g2} are considered. 
Finally, \code{g1} can be given as a matrix with two columns and each row indicating a pair. 

If \code{dat} is a list of \code{linkdat} objects, then all ID labels must be characters of the form "X-Y",
where X is the family ID, and Y the individual ID.
}
  \item{dat}{
Either a data frame (typically produced by \code{paramlink:::.ibdPrep}), a single \code{linkdat} object 
or a list of such.
}
  \item{f}{
A numeric vector containing the frequencies of the "2" allele. If NULL (default), this is computed 
from \code{dat}.
}
  \item{error}{
A single numeric to be used in error modelling, or NULL. See details.
}
  \item{plot.action}{
An integer: 2 = make new triangle and draw points; 1 = add points to existing triangle; 0 = no plotting. 
}
  \item{pointcol}{
A single numeric or character: Color specification for the points in the plot.
}
  \item{cex}{
A single numeric: Point size. (See\code{\link{par}})
}
  \item{\dots}{
Other arguments.
}
}
\details{
Maximum likelihood estimation is done by the \code{\link{maxLik}} function in the \code{maxLik} package.

More to come here.
}
\value{
A matrix (or a list of matrices) with 6 columns: ID1, ID2, N (#markers with no missing alleles), k0, k1, k2. 
}
\references{
E.A. Thompson (2000): Statistical Inferences from Genetic Data on Pedigrees. 
NSF-CBMS Regional Conference Series in Probability and Statistics. Volume 6.
}
\author{
Magnus Dehli Vigeland
}

\seealso{
\code{\link{maxLik}}, \code{\link{linkdat}}, \code{\link{points}}
}
\examples{
x = cousinPed(2)
x = simpleSim(x, 1000, 1:2) # 1000 equifrequent SNPs
examineKinships(x, who=c("par", "sib"))

\dontrun{
# comparing the founders with another unrelated person:
y = singleton(1, famid=2)
y = simpleSim(y, 1000, 1:2)
examineKinships(list(x, y), who="unrelated", interfam="founders")

# IBD estimates for the cousins (famid=1, iid=7 and 8) vs. 
# unrelated singleton (famid=2, iid=1):
ibd = IBDestimate(c("1-7", "1-8"), "2-1", list(x,y), error=0.005)
ibd
}
}
\keyword{math}

\name{markerSim}
\alias{markerSim}
\title{
Simulate SNP markers
}
\description{
Simulates SNP markers conditional on the pedigree structure, affection statuses and disease model.
}
\usage{
markerSim(x, N=1, available=x$orig.ids, alleles=NULL, afreq=NULL, partialmarker=NULL, 
          loop_breakers=NULL, eliminate=0, seed=NULL, method=3, verbose=TRUE)
}
\arguments{
  \item{x}{
a \code{linkdat} object
}
  \item{N}{
a positive integer: the number of markers to be simulated
}
 \item{available}{
a vector containing IDs of the available individuals, i.e. those whose genotypes should be simulated.
}
  \item{alleles}{
a vector containing the alleles for the marker to be simulation. If a single integer is given, this is interpreted as the number of alleles, and the actual alleles are 1:alleles. Should be NULL if \code{partialmarker} is not NULL.
}
  \item{afreq}{
a vector of length 2 containing the population frequencies for the marker alleles. Should be NULL if \code{partialmarker} is not NULL.
}
  \item{partialmarker}{
Either NULL (resulting in unconditional simulation), or a marker object on which the simulation should be conditioned. 
}
 \item{loop_breakers}{
a numeric containing IDs of individuals to be used as loop breakers. Relevant only if the pedigree has loops, and only if \code{partialmarker} is not NULL. See \code{\link{breakLoops}}.
}
  \item{eliminate}{
A non-negative integer, indicating the number of iterations in the internal genotype-compatibility algorithm. Positive values can save time if partialmarkers is non-NULL and the number of alleles is large.
}
  \item{seed}{
NULL, or a numeric seed for the random number generator. 
}
  \item{method}{
either 1,2 or 3: 
If 1, genotypes for all available individuals are simulated conditionally by brute force, one by one. This requires many likelihood calculations and is potentially very slow.
If 2, the joint genotype probabilities for some of the individuals are precomputed (optimized to ensure the minimum number of calls to the likelihood function). After sampling genotypes for these individuals, the rest are done by brute force as in method 1.
If 3, only the founders that are dependant of the partial marker data are simulated conditionally (using method 2), while the others are filled in last, by simple gene dropping. 
In rare cases this could involve brute force simulation of more individuals that the \code{available}, but should usually be the quickest.

If \code{partialmarker} is NULL, genotypes of all individuals are simulated by simple gene dropping, leaving the \code{method} option irrelevant.
}
  \item{verbose}{
a logical. 
}

}
\details{
This implementats (with various time savers) the algortihm used in SLINK of the LINKAGE/FASTLINK suite.
}
\value{
a \code{linkdat} object equal to \code{x} in all respects except its \code{markerdata} entry, which consists of the \code{N} simulated markers.
}
\author{
Magnus Dehli Vigeland
}
\seealso{
\code{\link{linkage.power}}
}
\references{
G. M. Lathrop, J.-M. Lalouel, C. Julier, and J. Ott, \emph{Strategies for Multilocus Analysis in Humans}, 
PNAS 81(1984), pp. 3443-3446.
}
\examples{
x = nuclearPed(2)
partial = marker(x, 3, 1, alleles=1:3)
markerSim(x, N=1, alleles=1:3)
markerSim(x, N=1, partialmarker=partial)
markerSim(x, N=1, partialmarker=partial, method=2)
markerSim(x, N=1, available=4, partialmarker=partial, method=3)
}
\keyword{ datagen }


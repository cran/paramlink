\name{simpleSim}
\alias{simpleSim}
\title{
Unconditional marker simulation
}
\description{
Unconditional simulation of unlinked markers
}
\usage{
simpleSim(x, N, alleles, afreq, available, Xchrom=FALSE, seed=NULL, verbose=T)
}
\arguments{
  \item{x}{
a \code{\link{linkdat}} object
}
  \item{N}{
a positive integer: the number of markers to be simulated
}
  \item{alleles}{
a vector containing the allele names. If missing, the alleles are taken to be \code{seq_along(afreq)}.
}
 \item{afreq}{
a vector of length 2 containing the population frequencies for the alleles. If missing, the alleles are assumed equifrequent.
}
 \item{available}{
a vector containing IDs of the available individuals, i.e. those whose genotypes should be simulated.
}
  \item{Xchrom}{
a logical: X linked markers or not?
}
  \item{seed}{
NULL, or a numeric seed for the random number generator. 
}
 \item{verbose}{
a logical.
}

}
\details{
This simulation is done by distributing alleles randomly to all founders, followed by unconditional gene dropping down throughout the pedigree (i.e. for each non-founder a random allele is selected from each of the parents). Finally the genotypes of any individuals not included in \code{available} are removed.
}
\value{
a \code{linkdat} object equal to \code{x} in all respects except its \code{markerdata} entry, which consists of the \code{N} simulated markers.
}
\author{
Magnus Dehli Vigeland
}
\seealso{
\code{\link{markerSim}}, \code{\link{linkageSim}}
}
\examples{
x = nuclearPed(1)
simpleSim(x, N=3, afreq=c(0.5, 0.5))

y = addOffspring(cousinPed(1), father=7, mother=8, noffs=1)
simpleSim(y, N=3, alleles=LETTERS[1:10])
}
\keyword{ datagen }
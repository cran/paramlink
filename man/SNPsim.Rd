\name{SNPsim}
\alias{SNPsim}
\title{
Simulate SNP markers
}
\description{
Simulates SNP markers conditional on the pedigree structure, affection statuses and disease model.
}
\usage{
SNPsim(x, N=1, partialmarker=NULL, available=x$available, afreq=c(0.5, 0.5), loop_breakers=NULL, 
       unique=FALSE, seed=NULL, verbose=TRUE)
}
\arguments{
  \item{x}{
a \code{linkdat} object
}
  \item{N}{
a positive integer: the number of SNPs to be simulated
}
  \item{partialmarker}{
Either NULL (indicating no given marker data), or a matrix with two columns giving alleles for each individual (0 = unkwown). 
}
 \item{available}{
a vector containing IDs of the available individuals, i.e. those whose genotypes should be simulated.
}
 \item{afreq}{
a vector of length 2 containing the population frequencies for the SNP alleles.
}
 \item{loop_breakers}{
a numeric containing IDs of individuals to be used as loop breakers. Relevant only if the pedigree has loops. See \code{\link{breakLoops}}.
}
  \item{unique}{
a logical indicating if duplicates among the simulated markers should be removed. 
}
  \item{seed}{
NULL, or a numeric seed for the random number generator. 
}
  \item{verbose}{
a logical. 
}

}
\details{
This is an implementation of the algortihm used in SLINK of the LINKAGE/FASTLINK suite.
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
data(toyped)
x = linkdat(toyped, model=1)
y = SNPsim(x, N=3)
plot(y, marker=1)

z = addOffspring(cousinPed(1), father=7, mother=8, noffs=1, aff=2)
z = setModel(z, 2)
SNPsim(z, N=1, loop_breakers=7)

m = marker(z, 9, c(2,2), alleles=1:2)
SNPsim(z, N=1, partialmarker=m, loop_breakers=7)

}
\keyword{ datagen }


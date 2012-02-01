\name{twoMarkerDistribution}
\alias{twoMarkerDistribution}
\title{
Genotype probability distribution
}
\description{
Computes the joint genotype distribution of two markers for a specified pedigree member, conditional on existing genotypes and pedigree information. 
Only for autosomal markers at the moment.
}
\usage{
twoMarkerDistribution(x, id, partialmarker1, partialmarker2, theta, 
                      loop_breakers=NULL, verbose=TRUE)
}

\arguments{
  \item{x}{
A \code{linkdat} object.
}
  \item{id}{
The individual in question.
}
  \item{partialmarker1, partialmarker2}{
Either a single integer indicating the number of one of \code{x}'s existing markers, or a \code{marker} object.
}
  \item{theta}{
A single numeric in the interval [0, 0.5] - the recombination fraction between the two markers.
}
  \item{loop_breakers}{
A numeric containing IDs of individuals to be used as loop breakers. Relevant only if the pedigree has loops. See \code{\link{breakLoops}}.
}
  \item{verbose}{
A logical.
}
}
\value{
A named matrix giving the joint genotype distribution.
}
\author{
Magnus Dehli Vigeland
}

\seealso{
\code{\link{oneMarkerDistribution}}
}
\examples{
data(toyped)
x = linkdat(toyped)
emptySNP = marker(x, 1:4, 0, alleles=c('a','b')) # empty marker
SNP1 = marker(x, 1, c(1,1), 2, c(1,0), alleles=1:2, afreq=c(0.1, 0.9))
twoMarkerDistribution(x, id=2, emptySNP, SNP1, theta=0)
twoMarkerDistribution(x, id=2, emptySNP, SNP1, theta=0.5)
twoMarkerDistribution(x, id=3, emptySNP, SNP1, theta=0.5)
}
\keyword{ distributions }
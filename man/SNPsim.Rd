\name{SNPsim}
\alias{SNPsim}
\title{
Simulate SNP markers
}
\description{
Simulates SNP markers conditional on the pedigree structure, affection statuses and disease model.
}
\usage{
SNPsim(x, N=1, unique=FALSE, seed=NULL)
}
\arguments{
  \item{x}{
a \code{linkdat} object
}
  \item{N}{
a positive integer: the number of SNPs to be simulated
}
  \item{unique}{
a logical indicating if duplicates among the simulated markers should be removed. 
}
  \item{seed}{
NULL, or a numeric seed for the random number generator. 
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
y = SNPsim(x, N=5)
plot(y, marker=1)
}
\keyword{ datagen }


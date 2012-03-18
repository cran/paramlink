\name{randomPed}
\alias{randomPed}
\title{
Random pedigree
}
\description{
Creates a random medical pedigree with specified number of generations. 
}
\usage{
randomPed(gen, lambda = 2, penetrances = c(0, 1, 1), 
          naff = "last.gen", founder.mut = 1)
}

\arguments{
  \item{gen}{
an integer in the interval \code{[2,5]} indicating the number of generations.
}
  \item{lambda}{
a positive numeric. For each descendent of the first generation, the number of offspring is sampled from a Poisson distribution with parameter \code{lambda}.
}
  \item{penetrances}{
a numeric of length 3, where for \code{i=1,2,3}, \code{penetrances[i] = Prob(affected | i-1 disease alleles)}.
}
  \item{naff}{
an integer specifying a lower bound on the number of affected individuals, or the character "last.gen". The latter produce a pedigree where at least one in the yongest generation is affected.
 }
  \item{founder.mut}{
an integer, the number of disease alleles to be distributed among the founders.
}
}
\details{
The function produce simple pedigrees. Each founder is given at most one disease allele. At least one of the two top founders carries a disease allele.
}
\value{
A \code{link{linkdat}} object.
}
\author{
Magnus Dehli Vigeland
}

\seealso{
\code{\link{linkdat}}
}
\examples{
plot(randomPed(3))

#gives error message: Not enough founder mutations
\dontrun{
randomPed(gen=4, penetrances=c(0,0,1), naff=2, founder.mut=1)
}
}
\keyword{ distributions }

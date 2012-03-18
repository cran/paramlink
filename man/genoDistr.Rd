\name{genoDistr}
\alias{genoDistr}
\title{
Genotype probability distribution
}
\description{
Computes genotype probabilities of a single individual, conditional on existing markers and pedigree information. 
Only for diallelic markers at the moment.
}
\usage{
genoDistr(x, id, partialmarker=NULL, t=0.5)
}

\arguments{
  \item{x}{
a \code{linkdat} object.
}
  \item{id}{
an integer indicating the individual in question.
}
  \item{partialmarker}{
either NULL (indicating no given marker data), a single integer indicating the number of an existing marker,
or a matrix with two columns giving alleles for each individual (0 = unkwown) 
}
  \item{t}{
a numeric, the recombination fraction between marker and disease locus
}
}
\value{
If \code{id} is male and the model is X-linked: A probability vector of length 3 (if the model is autosomal or \code{id} is female) or 2 (if the model is X-linked and \code{id} is male), giving the probabilities for the different genotypes.
}
\author{
Magnus Dehli Vigeland
}

\seealso{
\code{\link{likelihood}}
}
\examples{
data(toyped)
x = linkdat(toyped, model=1)
genoDistr(x, id = 1)

m1 = rbind(c(1,1),c(1,0),c(0,0),c(0,0))
genoDistr(x, id = 3, partialmarker = m1)

m2 = rbind(c(0,0),c(0,0),c("A","B"),c("A","B"))
genoDistr(x, id = 1, partialmarker = m2)

\dontrun{
# gives Mendelian error (correctly):
m3 = rbind(c(1,1),c(1,1),c(2,2),c(0,0))
genoDistr(x, id = 4, partialmarker = m3)
}

# X-linked example
data(Xped)
y = linkdat(Xped, model=4) #X-linked recessive model
# create empty marker, and then give individuals 5 and 15 genotype AA:
y = setMarkers(y, matrix(0, ncol=2, nrow=15))
y = modifyMarker(y, id=c(5, 15), alleles=c('A', 'A'))	
genoDistr(y, id=13, partialmarker=1)	#results: A - 0.8; B - 0.2
genoDistr(y, id=13, partialmarker=1, t=0)	#results: A - 0.667; B - 0.333
}

\keyword{ distributions }
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
genoDistr(x, id, partialmarker=NULL, t=0.5, original.id=TRUE)
}

\arguments{
  \item{x}{
A \code{linkdat} object.
}
  \item{id}{
The individual in question.
}
  \item{partialmarker}{
Either NULL (indicating no given marker data), a single integer indicating the number of an existing marker,
or a matrix with two columns giving alleles for each individual (0 = unkwown). 
}
  \item{t}{
The recombination fraction between marker and disease locus.
}
  \item{original.id}{
A logical indicating whether 'id' refers to the original IDs or the internal labeling.
}
}
\value{
A numeric vector giving the probabilities for the different genotypes. The vector has length 3 (corresponding to the genotypes, e.g. AA, AB and BB) except if the model is X-linked and \code{id} is male. In that case the output has length 2 (corresponding to the two alleles A and B).
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
genoDistr(x, id=1)

m1 = rbind(c(1,1),c(1,0),c(0,0),c(0,0))
genoDistr(x, id=3, partialmarker=m1)

m2 = rbind(c(0,0),c(0,0),c("A","B"),c("A","B"))
genoDistr(x, id=1, partialmarker=m2)

\dontrun{
# gives Mendelian error (correctly):
m3 = rbind(c(1,1),c(1,1),c(2,2),c(0,0))
genoDistr(x, id=4, partialmarker=m3)}

# another example, the genotype distribution of an individual whose half cousin is homozygous for a rare allele.
y = halfCousinPed(degree=1)
y = setModel(y, 1, afreq=c(0.99, 0.01)) #autosomal model, allele frequencies for the marker set to 99% and 1%.
y = setMarkers(y, cbind(c(0,0,0,0,0,0,0,0,'b'), c(0,0,0,0,0,0,0,0,'b')))
plot(y, marker=1)
genoDistr(y, id=8, partialmarker=1)

# X-linked example
data(Xped)
z = linkdat(Xped, model=4) #X-linked recessive model
# create empty marker, and then give individuals 5 and 15 genotype AA:
z = setMarkers(z, matrix(0, ncol=2, nrow=15))
z = modifyMarker(z, id=c(5, 15), alleles=c('A', 'A'))	
genoDistr(z, id=13, partialmarker=1)	#results: A - 0.8; B - 0.2
genoDistr(z, id=13, partialmarker=1, t=0)	#results: A - 0.667; B - 0.333
}

\keyword{ distributions }
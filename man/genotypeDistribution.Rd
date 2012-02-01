\name{oneMarkerDistribution}
\alias{oneMarkerDistribution}
\title{
Genotype probability distribution
}
\description{
Computes the (joint) genotype probability distribution of one or several pedigree members, possibly conditional on partial marker data.
Only for autosomal markers at the moment.
}
\usage{
oneMarkerDistribution(x, ids, partialmarker, theta=0.5, loop_breakers=NULL, 
                     verbose=TRUE)
}

\arguments{
  \item{x}{
A \code{linkdat} object.
}
  \item{ids}{
A numeric with ID labels of one or more pedigree members.
}
  \item{partialmarker}{
Either a single integer indicating the number of one of \code{x}'s existing markers, or a \code{\link{marker}} object.
}
  \item{theta}{
The recombination fraction between marker and disease locus. If \code{theta < 0.5} an error is given if \code{x} does not have a disease model (see \code{\link{setModel}}).
If \code{theta = 0.5} (the default), the marker is assumed to segregate independently of the disease, and there is no need to specify a disease model for \code{x}. 
}
  \item{loop_breakers}{
A numeric containing IDs of individuals to be used as loop breakers. Relevant only if the pedigree has loops. See \code{\link{breakLoops}}.
}
  \item{verbose}{
A logical.
}
}
\value{
A named array (of dimension \code{length(ids)}) giving the joint genotype distribution for the \code{ids} individuals, conditional on the non-missing alleles in \code{partialmarker}.
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
	
snp1 = marker(x, 1, c(1,1), 2, c(1,0), alleles=1:2, afreq=c(0.1, 0.9))
oneMarkerDistribution(x, ids=3, partialmarker=snp1)
oneMarkerDistribution(x, ids=3:4, partialmarker=snp1)


snp2 = marker(x, 3:4, c("A","B"))
oneMarkerDistribution(x, ids=1, partialmarker=snp2)
oneMarkerDistribution(x, ids=1, partialmarker=snp2, theta=0)

\dontrun{
# gives Mendelian error (correctly):
snp.err = marker(x, allelematrix=rbind(c(1,1),c(1,1),c(2,2),c(0,0)))
oneMarkerDistribution(x, ids=4, partialmarker=snp.err)
mendelianCheck(setMarkers(x, snp.err))
}

# another example, the genotype distribution of an individual whose half cousin is homozygous for a rare allele.
y = halfCousinPed(degree=1)
snp = marker(y, 9, c('a','a'), alleles=c('a', 'b'), afreq=c(0.01, 0.99))
oneMarkerDistribution(y, ids=8, partialmarker=snp)

\dontrun{
# X-linked example, doesn't work in this version.
data(Xped)
z = linkdat(Xped, model=4) # X-linked recessive model
snpX = marker(z, c(5,15), c('A','A'), alleles=c('A', 'B'))	
res_0.5 = oneMarkerDistribution(z, ids=13, partialmarker=snpX, theta=0.5)	# results: A - 0.8; B - 0.2
res_0.0 = oneMarkerDistribution(z, ids=13, partialmarker=snpX, theta=0)		# results: A - 0.667; B - 0.333
res_0.5
res_0.0
stopifnot(round(res_0.5[1], 2)==0.8, round(res_0.0[1], 2) == 0.67)
}
}
\keyword{ distributions }
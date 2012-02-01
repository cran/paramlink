\name{likelihood}
\alias{likelihood}
\title{Pedigree likelihood
}
\description{
Calculates the likelihood of a pedigree, given the recombination rate between the disease locus and one marker locus.
}
\usage{
likelihood(x, marker, theta=NULL, afreq=NULL, logbase=NULL, TR.MATR=NULL, 
           singleNum.geno=NULL)
}
\arguments{
  \item{x}{
a \code{linkdat} object, with \code{x$model} different from NULL.
}
  \item{marker}{
an integer between 0 and \code{x$nMark}, indicating which marker to use in the calculation. If NULL, the marker is taken to be all zeroes.
}
  \item{theta}{
a numeric of length 1 indicating the recombination rate between the disease and marker loci. 
To make biological sense this should be between 0 and 0.5.
}
  \item{afreq}{
a numeric containing the marker allele frequencies. 
}
 \item{logbase}{
a numeric, or NULL. If numeric, the loglikelihood is computed, with \code{logbase} as basis for the logarithm.
}
 \item{TR.MATR}{
an 10*10*10 array with haplotype transmission probabilities, for speedy computation with diallelic markers.
	}
 \item{singleNum.geno}{
Mostly for internal use: a numeric vector containing genotypes coded as a single integer: 00 -> 0, 11 -> 1, 22 -> 2, 12/21 -> 3, 01/10 -> 4, 02/20 -> 5. 
If this is non-null the \code{marker} argument is ignored.
}
}
\details{
The likelihood is calculated using the Elston-Stewart algorithm.
}
\value{
The likelihood of the pedigree with disease marker data, given the disease model and the recombination rate. 
If \code{marker=NULL}, the result is the likelihood of the pedigree at the disease locus only (regardless the value of \code{theta}).

If the parameter \code{logbase} is numeric, the output is \code{log(likelihood, logbase)}.
}
\author{
Magnus Dehli Vigeland
}

\seealso{
\code{\link{lod}}
}
\examples{
ped = data.frame(ID=1:4, FID=c(0,0,1,1), MID=c(0,0,2,2), SEX=c(1,2,1,1), 
                 AFF=c(2,1,2,2), M_A1=c(1,1,1,1), M_A2=c(2,1,2,2))

x = linkdat(ped, model = 1)


lod1 = likelihood(x, marker=1, theta=0, logbase=10) - 
       likelihood(x, marker=1, theta=0.5, logbase=10)
lod2 = lod(x, markers=1, theta=0)
# should be the same:
stopifnot(identical(lod1, as.numeric(lod2)))
}

\keyword{ math }

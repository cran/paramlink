\name{likelihood}
\alias{likelihood}
\alias{likelihood.linkdat}
\alias{likelihood.singleton}
\alias{likelihood.list}
\alias{likelihoodSNP}
\title{Pedigree likelihood
}
\description{
Calculates various forms of pedigree likelihoods.
}
\usage{
likelihood(x, ...)
\method{likelihood}{linkdat}(x, locus1, locus2=NULL, theta=NULL, 
                    startdata=NULL, eliminate=0, logbase=NULL, ...)
\method{likelihood}{singleton}(x, locus1, logbase=NULL, ...)
\method{likelihood}{list}(x, locus1, ...)

likelihoodSNP(x, marker, theta=NULL, afreq=NULL, logbase=NULL, TR.MATR=NULL, 
              singleNum.geno=NULL)
}
\arguments{
  \item{x}{
a \code{linkdat} object, a \code{singleton} object, or a list of such objects. In \code{likelihoodSNP}, 
\code{x} must be a \code{linkdat} object, with \code{x$model} different from NULL.
}
  \item{locus1}{
a \code{\link{marker}} object compatible with \code{x}. If \code{x} is a list, then \code{locus1} must be a list of corresponding
\code{marker} objects.
}
  \item{locus2}{
either NULL, the character "disease", or a \code{\link{marker}} object compatible with \code{x}. See Details.
}
  \item{theta}{
the recombination rate between locus1 and locus2 (in \code{likelihoodSNP}: between the marker and the disease locus). 
To make biological sense theta should be between 0 and 0.5.
}
  \item{eliminate}{
mostly for internal use: a non-negative integer indicating the number of iterations in the internal genotype-compatibility algorithm. 
Positive values can save time if \code{partialmarker} is non-empty and the number of alleles is large.
}
  \item{logbase}{
a numeric, or NULL. If numeric the loglikelihood is returned, with \code{logbase} as basis for the logarithm.
}
  \item{marker}{
an integer between 0 and \code{x$nMark}, indicating which marker to use in the calculation.
}
  \item{afreq}{
a numeric containing the marker allele frequencies. 
}
  \item{startdata}{
for internal use.
}
 \item{TR.MATR}{
for internal use: an 10*10*10 array with haplotype transmission probabilities, for speedy computation with diallelic markers.
}
  \item{singleNum.geno}{
for internal use: a numeric vector containing SNP genotypes coded as a single integer: 00 -> 0, 11 -> 1, 22 -> 2, 12/21 -> 3, 01/10 -> 4, 02/20 -> 5. 
If this is non-NULL the \code{marker} argument is ignored.
}
  \item{\dots}{
further arguments. Not used.
}
}
\details{
The likelihoods are calculated using the Elston-Stewart algorithm.
If \code{locus2 = NULL}, the result is the likelihood of the genotypes observed at the marker in locus1.
If \code{locus2 = "disease"}, the result is the likelihood of the marker genotypes in locus1, given the affection statuses of the pedigree members, the disease model and the recombination rate \code{theta} between the marker and disease loci. 
(This computation is central in parametric linkage analysis.)
If \code{locus2} is a marker object, the result is the likelihood of the genotypes at the two markers, given the recombination rate theta between them. 

The function \code{likelihoodSNP} is a fast version of \code{likelihood.linkdat} in the case where \code{locus2 = "disease"} and the marker in locus1 is diallelic.
}
\value{
The likelihood of the data. If the parameter \code{logbase} is numeric, the output is \code{log(likelihood, logbase)}.
}
\author{
Magnus Dehli Vigeland
}

\seealso{
\code{\link{lod}}
}
\examples{
data(toyped)
x = linkdat(toyped, model=1) #dominant model

lod1 = likelihoodSNP(x, marker=1, theta=0, logbase=10) - 
       likelihoodSNP(x, marker=1, theta=0.5, logbase=10)
lod2 = lod(x, markers=1, theta=0)
# should be the same:
stopifnot(identical(lod1, as.numeric(lod2)), round(lod1, 2)==0.3)
}

\keyword{ math }

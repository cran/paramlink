\name{setModel}
\alias{setModel}
\alias{print.linkdat.model}

\title{
Set, change or display the model parameters for 'linkdat' objects
}
\description{
Functions to set, change and display the model parameters involved in parametric linkage analysis.
}
\usage{
setModel(x, model = NULL, chrom = NULL, penetrances = NULL, 
         dfreq = NULL, nallel = NULL, afreq = NULL)

\method{print}{linkdat.model}(x, ...)
}

\arguments{
  \item{x}{
in \code{setModel}: a \code{linkdat} object. In \code{print.linkdat.model}: a \code{linkdat.model} object, as described below. 
}
  \item{model}{
NULL, or an object of class \code{linkdat.model}, namely a list with elements 
\code{chrom}, \code{penetrances}, \code{dfreq}, \code{nallel}, \code{afreq}. 
In the \code{setModel} function, the \code{model} parameter can also be one of the integers 1-4, with the following meanings:

1 = autosomal dominant, i.e. \code{list(chrom="AUTOSOMAL", penetrances=c(0,1,1), dfreq=1e-5, nallel=2, afreq=c(.5,.5))}

2 = autosomal recessive, i.e. \code{list(chrom="AUTOSOMAL", penetrances=c(0,0,1), dfreq=1e-5, nallel=2, afreq=c(.5,.5))}

3 = X-linked dominant, i.e. \code{list(chrom="X", penetrances=list(male=c(0,1), female=c(0,1,1)), dfreq=1e-5, nallel=2, afreq=c(.5,.5))}

4 = X-linked recessive, i.e. \code{list(chrom="X", penetrances=list(male=c(0,1), female=c(0,0,1)), dfreq=1e-5, nallel=2, afreq=c(.5,.5))}
}
  \item{chrom}{
a character, either 'AUTOSOMAL' or 'X'. Lower case versions are allowed and will be converted automatically.
}
  \item{penetrances}{
if \code{chrom=='AUTOSOMAL'}: a numeric of length 3 - \code{(f_0, f_1, f_2)} - where \code{f_i} is the probability of being affected given \code{i} disease alleles. 
If \code{chrom=='X'}: a list of length two containing the penetrances for each sex: \code{penetrances = list(male=c(f0, f1), female=c(f0, f1, f2))}. 
}
  \item{dfreq}{
the population frequency of the disease allele.
}
  \item{nallel}{
the number of alleles at the marker loci. At the moment only diallelic SNPs are allowed in likelihood calculations.
}
  \item{afreq}{
a numeric of length \code{nallel} - the population allele frequencies at the marker loci. 
}
  \item{...}{further parameters}
}

\value{
a \code{linkdat} object equal to \code{x}, except that its \code{model} entry is set according to the given parameters. 
}

\author{
Magnus Dehli Vigeland
}
\seealso{
\code{\link{linkdat}}
}
\examples{
data(toyped)
x = linkdat(toyped)
x = setModel(x, model=1)
summary(x)
y = setModel(x, model=4, dfreq=0.1)
summary(y)
}
\keyword{ math }

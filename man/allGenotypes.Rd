\name{allGenotypes}
\alias{allGenotypes}
\alias{fast.grid}
\alias{geno.grid.subset}
\title{
Genotype combinations
}
\description{
Auxiliary functions computing possible genotype combinations in a pedigree. These are not normally intended for end users.
}
\usage{
geno.grid.subset(x, partialmarker, ids, chrom, make.grid=T)

allGenotypes(n)

fast.grid(argslist, as.list=FALSE)
}
\arguments{
  \item{x}{
a \code{\link{linkdat}} object.
}
  \item{partialmarker}{
a \code{\link{marker}} object compatible with \code{x}.
}
  \item{ids}{
a numeric with ID labels of one or more pedigree members.
}
  \item{chrom}{
a character, either 'X' or 'AUTOSOMAL'. If missing, the 'chrom' attribute of \code{partialmarker} is used. If this is also missing, then 'AUTOSOMAL' is taken as the default value.
}
  \item{make.grid}{
a logical. If FALSE, a list is returned, otherwise \code{fast.grid} is applied to the list before returning it.
}
  \item{n}{
a positive integer.
}
  \item{argslist}{
a list of vectors.
}
  \item{as.list}{
if TRUE, the output is a list, otherwise a matrix.
}
}
\value{
\code{allGenotypes} returns a matrix with 2 columns and \code{n + n*n(n-1)/2} rows containing all possible (unordered) genotypes at a biallelic locus with alleles \code{1,2,\dots,n}.
\code{fast.grid} is basically a stripped down version of \code{\link{expand.grid}}.
}
\author{
Magnus Dehli Vigeland
}
\examples{
m = allGenotypes(2)
stopifnot(m == rbind(c(1,1), c(2,2), 1:2))
}
\keyword{ methods }
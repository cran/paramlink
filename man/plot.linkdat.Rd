\name{plot.linkdat}
\alias{plot.linkdat}
\title{
Plot pedigrees with genotypes
}
\description{
Plots a medical pedigree with genotypes for a single marker, or (for power analysis) indicating which individuals are available for genotyping.  
}
\usage{
\method{plot}{linkdat}(x, marker = NULL, alleles = NULL, sep = "", missing = "-", 
                       id.labels = x$orig.ids, title = paste("Family", x$famid), 
                       available = FALSE, col =1, deceased = numeric(0), 
                       starred = numeric(0), aff2 = NULL, margins = c(4.1, 1, 4.1, 1), ...)
}

\arguments{
  \item{x}{
a \code{linkdat} object.
}
  \item{marker}{
NULL, or a single integer.
}
  \item{alleles}{
a character vector with allele names.
}
  \item{sep}{
a character of length 1 separating alleles for diploid markers. 
}
  \item{missing}{
the symbol (integer or character) for missing alleles.
}
  \item{id.labels}{
a vector with labels/names of the pedigree members.
}
  \item{title}{
the plot title. If NULL or "", no title is added to the plot. 
}
  \item{available}{
if TRUE, the available individuals are drawn in red. 
}
  \item{col}{
a vector with color indicators for the pedigree members. Recycled if necessary. By default everyone is drawn black.
}
  \item{deceased}{
a numeric containing ID's of deceased pedigree members. 
}
  \item{starred}{
a numeric containing ID's of pedigree members that should be marked with a star in the pedigree plot.
}
  \item{aff2}{
NULL, or a numeric with affection statuses (2=affected, 1=unaffected, 0=unknown) for a second trait.
}
  \item{margins}{
a numeric of length 4 indicating the plot margins.
}
  \item{\dots}{
arguments to \code{plot.pedigree} from the \code{kinship} package. In particular \code{symbolsize} is useful.}
}
\details{
Uses the \code{plot.pedigree} from the \code{kinship} package.
}

\author{
Magnus Dehli Vigeland
}

\seealso{
\code{\link{plot.pedigree}}}
\examples{
data(toyped)
x = linkdat(toyped)
plot(x, marker=1, alleles=c("a1","a2"), sep=" | ", deceased=2)
}

\keyword{ plot }

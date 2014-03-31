\name{as.data.frame.linkdat}
\alias{as.data.frame.linkdat}
\title{
linkdat to data.frame conversion
}
\description{
Convert a linkdat object to data.frame for pretty printing.
}
\usage{
\method{as.data.frame}{linkdat}(x, ..., famid=F, markers=seq_len(x$nMark), alleles=NULL, 
              missing=NULL, singleCol=FALSE, sep="")
}
\arguments{
 \item{x}{
a \code{\link{linkdat}} object.
}
 \item{famid}{
a logical indicating if the family identifier should be included as the first column.
}
 \item{markers}{
a numeric indicating which markers should be included/printed.
}
 \item{alleles}{
a character containing allele names, e.g. \code{alleles=c("A","B")}.
}
 \item{missing}{
the character (of length 1) used for missing alleles. Defaults to "0".
}
 \item{singleCol}{
a logical: Should the two alleles for each marker be pasted into one column or kept in separate columns?
}
 \item{sep}{
a single character to be used as allele separator if \code{singleCol=TRUE}.
}
 \item{\dots}{
further arguments (not used).
}
} 
\details{
This function is mainly intended for pretty-printing \code{linkdat} objects (for instance it is called by \code{\link{print.linkdat}}).
For direct manipulation of the pedigree and/or marker matrices, it is better to use \code{\link{as.matrix.linkdat}}.
}
\value{
A \code{data.frame}.
}

\author{
Magnus Dehli Vigeland
}

\seealso{
\code{\link{print.linkdat}}, \code{\link{as.matrix.linkdat}}
}
\examples{
x = linkdat(toyped)
x

# Printing x as above is equivalent to:
as.data.frame(x, sep = "/", missing = "-", singleCol = TRUE)
}

\keyword{ IO }
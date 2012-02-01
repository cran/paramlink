\name{linkres}
\alias{linkres}
\alias{print.linkres}
\alias{summary.linkres}
\alias{plot.linkres}
\title{
S3 methods for class 'linkres'.
}
\description{
Functions for printing, summarizing and plotting the results of a linkage analysis.
}
\usage{
\method{print}{linkres}(x, ...)

\method{summary}{linkres}(object, threshold=NULL, ...)

\method{plot}{linkres}(x, chrom=NULL, ylim = NULL, ...)
}
\arguments{
  \item{x, object}{
a \code{linkres} object (normally produced by \code{\link{lod}} or \code{\link{merlin}}).
}
  \item{threshold}{
NULL, or a positive real number. In the latter case the summary will include the percentage of simulated markers whose LOD score exceeds this number.
}
  \item{chrom}{
NULL, or a numeric containing chromosome numbers. In the latter case only results for the markers on the indicated chromosomes will be plotted.
}
 \item{ylim}{
NULL, or a numeric of length 2: to be passed on to plot.default.
}
  \item{\dots}{
further arguments.
}
}

\author{
Magnus Dehli Vigeland
}

\seealso{
\code{\link{lod}}, \code{\link{merlin}}
}
\examples{
data(toyped)
x = linkdat(toyped, model=1)
lods = lod(x, theta="max")
summary(lods)
}
\keyword{ methods }

\name{linkres}
\alias{linkres}
\alias{as.data.frame.linkres}
\alias{print.linkres}
\alias{summary.linkres}
\alias{plot.linkres}
\alias{peakSummary}
\title{
S3 methods for class 'linkres'.
}
\description{
Functions for printing, summarizing and plotting the results of a linkage analysis.
}
\usage{
\method{as.data.frame}{linkres}(x, ..., sort=TRUE)

\method{plot}{linkres}(x, chrom=NULL, ylim = NULL, ...)

\method{print}{linkres}(x, ...)

\method{summary}{linkres}(object, ...)

peakSummary(x, threshold, width=1, physmap=NULL)
}
\arguments{
  \item{x, object}{
a \code{linkres} object (normally produced by \code{\link{lod}} or \code{\link{merlin}}).
}
  \item{sort}{
a logical, indicating if the data frame should be sorted according to map position.
}
  \item{chrom}{
NULL, or a numeric containing chromosome numbers. In the latter case only results for the markers on the indicated chromosomes will be plotted.
}
  \item{ylim}{
NULL, or a numeric of length 2: to be passed on to plot.default.
}
 \item{threshold}{
a single numeric. A peak is defined as a regions of at least \code{width} consecutive markers LOD score above \code{threshold}.
} 
  \item{width}{
a single numeric.
} 
  \item{physmap}{
a matrix or data frame with three columns: Marker name, chromosome and physical position. This argument is optional.
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
x = linkdat(toyped, model=1)
lods = lod(x, theta="max")
summary(lods)
as.data.frame(lods)
}
\keyword{ methods }
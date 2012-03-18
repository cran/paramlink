\name{linkres}
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

\method{plot}{linkres}(x, ylim = NULL, ...)
}
\arguments{
  \item{x, object}{
a \code{linkres} object (normally produced by \code{\link{lod}} or \code{\link{linkage.power}}.
}
	\item{threshold}{
a positive real number, or NULL. If non-null, the summary will include the percentage of simulated markers whose LOD score exceeds \code{threshold}.
	}
 \item{ylim}{
NULL, or a numeric of length 2; to be passed on to plot.default.
}
  \item{\dots}{
further arguments.
}
}

\author{
Magnus Dehli Vigeland
}

\seealso{
\code{\link{lod}},\code{\link{linkage.power}}
}
\examples{
data(toyped)
x=linkdat(toyped, model=1)
lods=lod(x, t="max")
summary(lods)
}
\keyword{ methods }

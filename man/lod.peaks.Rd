\name{lod.peaks}
\alias{lod.peaks}
\title{
LOD score peaks
}
\description{
Identify LOD score peaks
}
\usage{
lod.peaks(x, threshold, width = 1)
}
\arguments{
  \item{x}{
a \code{\link{linkres}} object
}
  \item{threshold}{
a single numeric
}
  \item{width}{
a positive integer
}
}
\details{
The function first transforms \code{x} to a data frame (using \code{\link{as.data.frame.linkres}} with \code{sort=T}. A peak is defined a run of at least \code{width}
consecutive markers with LOD score above or equal to \code{threshold}. If possible, one flanking marker is included on each side of the peak.
}
\value{
A list of data frames.
}
\author{
MAgnus Dehli Vigeland
}
\seealso{
\code{\link{linkres}}, \code{\link{lod}}, \code{\link{merlin}},
}
\examples{
## minimal example
data(toyped)
x = linkdat(toyped, model=1)
res = lod(x)
peak1 = lod.peaks(res, threshold=0)
peak2 = lod.peaks(res, threshold=0, width=2)
peak3 = lod.peaks(res, threshold=1)
stopifnot(length(peak1)==1, nrow(peak1[[1]])==1, length(peak2)==0, length(peak3)==0)
}
\keyword{ manip }


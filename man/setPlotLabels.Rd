\name{setPlotLabels}
\alias{setPlotLabels}
\title{
Attach plot labels to a linkdat object
}
\description{
This function attaches (or modifies) a character vector of plotting labels for the pedigree members of a linkdat object. 
This is useful since only numerical ID's are allowed in defining pedigrees in paramlink.
}
\usage{
setPlotLabels(x, labels, ids = x$orig.ids)
}
\arguments{
  \item{x}{
A linkdat object.
}
  \item{labels}{
A character vector of the same length as \code{ids}.
}
  \item{ids}{
A numeric vector of numerical IDs. Must be a subset of \code{x$orig.ids}.
}
}
\value{
A new linkdat object, differing from \code{x} only in \code{x$plot.labels}.
}
\author{
Magnus Dehli Vigeland
}
\seealso{
\code{\link{plot.linkdat}}
}
\examples{
x =  nuclearPed(1)
x = setPlotLabels(x, labels=c("Father", "Mother", "Son"))
plot(x)
}
\keyword{ dplot }

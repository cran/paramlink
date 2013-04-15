\name{setAvailable}
\alias{setAvailable}
\alias{swapAvailable}
\title{Functions for modifying availability vectors}
\description{
Functions to set and modify the availability vector of a 'linkdat' object. This vector is used in 'linkage.power' and 'linkageSim', indicating for whom genotypes should be simulated. 
}
\usage{
setAvailable(x, available)

swapAvailable(x, ids)
}
\arguments{
  \item{x}{
a \code{\link{linkdat}} object
}
 \item{available}{
a numeric containing the IDs of available individuals.
} 
 \item{ids}{
the individual(s) whose availability status should be swapped.
}
}
\value{
The modified \code{linkdat} object.
}
\author{
Magnus Dehli Vigeland
}
\seealso{
\code{\link{plot.linkdat}}, \code{\link{linkage.power}}, \code{\link{linkageSim}}
}
\examples{
data(toyped)
x = linkdat(toyped)
x = setAvailable(x, 3:4)
x = swapAvailable(x, 2:3)
x$available
}
\keyword{ manip }

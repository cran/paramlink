\name{setSim}
\alias{setSim}
\alias{swapSim}
\title{Simulation vector for 'linkdat' objects}
\description{
Functions to set and modify the simulation vector of a 'linkdat' object. This vector is used in 'linkage.power' and 'SNPsim', indicating for whom genotypes should be simulated. 
}
\usage{
setSim(x, simstatus)

swapSim(x, ids)
}
\arguments{
  \item{x}{
a \code{\link{linkdat}} object
}
 \item{simstatus}{
either NULL, a logical, or a numeric consisting of 0's and 2's (where 2=included and 0=excluded in the simulation). Recycled if necessary.
} 
 \item{ids}{
the individual(s) whose simulation status should be swapped.
}
}
\details{
If the simulation vector of a linkdat object \code{x} has been set, the included individuals are marked by a star in the pedigree plot.
}
\value{
The modified \code{linkdat} object.
}
\author{
Magnus Dehli Vigeland
}
\seealso{
\code{\link{linkage.power}}, \code{\link{SNPsim}}
}
\examples{
data(toyped)
x = linkdat(toyped)
x = setSim(x, TRUE)
x = swapSim(x, 1:2)
}
\keyword{ manip }

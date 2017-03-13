\name{is.linkdat}
\alias{is.linkdat}
\alias{is.singleton}
\alias{is.linkdat.list}
\title{
Is an object a linkdat object?
}
\description{
Functions for checking whether an object is a linkdat object, a singleton or a list of such.
}
\usage{
is.linkdat(x)
is.singleton(x)
is.linkdat.list(x)
}
\arguments{
  \item{x}{
Any \code{R} object.
}
}
\details{
Note that the \code{singleton} class inherits from \code{linkdat}, so if \code{x} is a singleton, \code{is.linkdat(x)} returns TRUE.
}
\value{
For \code{is.linkdat}: TRUE if \code{x} is a linkdat (or singleton) object, and FALSE otherwise.\cr
For \code{is.singleton}: TRUE if \code{x} is a singleton object, and FALSE otherwise.\cr
For \code{is.linkdat.list}: TRUE if \code{x} is a list of linkdat/singleton objects.
}
\author{
Magnus Dehli Vigeland
}
\seealso{
\code{\link{linkdat}}
}
\examples{
x1 = nuclearPed(1)
x2 = singleton(1)
stopifnot(is.linkdat(x1), !is.singleton(x1), 
          is.linkdat(x2), is.singleton(x2),
          is.linkdat.list(list(x1,x2)))
}
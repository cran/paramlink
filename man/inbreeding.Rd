\name{inbreeding}
\alias{inbreeding}
\alias{jacquard}
\title{
Inbreeding coefficients and Jacquard coefficients
}
\description{
Computes inbreeding coefficients for all pedigree members, and Jaquard's condensed identity coefficients for any pair of members. These are simple wrappers for functions in other packages, respectively \code{kinship2} and \code{identity}.
}
\usage{
inbreeding(x)
jacquard(x, ids)
}
\arguments{
  \item{x}{
a \code{\link{linkdat}} object
}
  \item{ids}{
a integer vector of length 2
}
}
\details{
The main work of \code{inbreeding} is done by \code{kinship2::kinship}. Similarly, \code{jacquard} wraps \code{identity::identity.coefs}.
}
\value{
For \code{inbreeding}, a numerical vector with the inbreeding coefficients, in the same order as \code{x$orig.ids}.
For \code{jaquard}, a numerical vector of length 9 (in the standard order of Jacquard's condensed coefficients).
}
\author{
Magnus Dehli Vigeland
}

\seealso{
\code{\link{kinship}},\code{identity::identity.coefs}
}
\examples{
# Offspring of first cousins
x = cousinsPed(1, child=TRUE)
inb = inbreeding(x)
stopifnot(inb[9] == 1/16)
}
\keyword{ math }

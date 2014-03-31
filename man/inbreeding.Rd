\name{inbreeding}
\alias{inbreeding}
\title{
Inbreeding coefficients
}
\description{
Computes inbreeding coefficients for all pedigree members.
}
\usage{
inbreeding(x)
}
\arguments{
  \item{x}{
a \code{\link{linkdat}} object
}
}
\details{
The main work is done by the \code{\link{kinship}} function the \code{kinship2} package.
}
\value{
A numerical vector with the inbreeding coefficients, in the same order as \code{x$orig.ids}.
}
\author{
Magnus Dehli Vigeland
}

\seealso{
\code{\link{kinship}},
}
\examples{
# Offspring of first cousins
x = addOffspring(cousinPed(1), father=7, mother=8, noffs=1)
inb = inbreeding(x)
stopifnot(inb[9] == 1/16)
}
\keyword{ math }

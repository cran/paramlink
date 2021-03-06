% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/relatednessCoeff.R
\name{relatednessCoeff}
\alias{relatednessCoeff}
\alias{inbreeding}
\alias{kinship.coeff}
\alias{jacquard}
\alias{jacquard2}
\title{Relatedness coefficients}
\usage{
inbreeding(x)

\method{kinship}{coeff}(x, ids = NULL)

jacquard(x, ids)

jacquard2(x, ids, verbose = FALSE, cleanup = TRUE)
}
\arguments{
\item{x}{a \code{\link{linkdat}} object.}

\item{ids}{a integer vector of length 2.}

\item{verbose}{a logical, indicating if messages from IdCoefs should be
printed.}

\item{cleanup}{a logical: If TRUE, the pedfile and sample file created for
the IdCoefs run are deleted automatically.}
}
\value{
For \code{inbreeding}, a numerical vector with the inbreeding
coefficients, with names according to the ID labels \code{x$orig.ids}.\cr
For \code{kinship.coeff}, either a single numeric (if \code{ids} is a pair
of pedigree members) or the whole kinship matrix, with \code{x$orig.ids} as
dimnames.\cr For \code{jaquard} and \code{jaquard2}, a numerical vector of
length 9 (in the standard order of Jacquard's condensed identity
coefficients).
}
\description{
Computes inbreeding coefficients for all pedigree members, and Jaquard's
condensed identity coefficients for any pair of members. These are simple
wrappers for functions in other packages or external programs.
}
\details{
Both \code{inbreeding} and \code{kinship.coeff} are thin wrappers of
\code{\link[kinship2]{kinship}}.  Similarly, \code{jacquard} wraps
\code{identity::identity.coefs}, which is an R implementation of the C
program \code{IdCoefs} written by Mark Abney. The function
\code{identity::identity.coefs} sometimes causes R to crash, hence as an
alternative wrapper, \code{jacquard2} executes an external call to the
original C program \code{IdCoefs} (version 2.1.1). Of course, this requires
\code{IdCoefs} to be installed on the computer (see link in the References
section below) and the executable to be in a folder included in the PATH
variable. The \code{jacquard2} wrapper works by writing the necessary files
to disk and calling \code{IdCoefs} via \code{\link{system}}.
}
\examples{

# Offspring of first cousins
x = cousinsPed(1, child=TRUE)
inb = inbreeding(x)
stopifnot(inb[9] == 1/16)

# if ID labels are not 1:9, care must be taken in extracting correct elements.
set.seed(1357)
y = relabel(x, sample(1:9))
child = leaves(y)
inbreeding(y)[child] #wrong
inb = inbreeding(y)[as.character(child)] #correct
inb
# the inbreeding coeff of the child equals the kinship coeff of parents
kin = kinship.coeff(y, parents(y, child))
stopifnot(inb==kin, inb==1/16)

}
\references{
The \code{IdCoefs} program: Abney, Mark (2009). A graphical
algorithm for fast computation of identity coefficients and generalized
kinship coefficients. Bioinformatics, 25, 1561-1563.
\url{http://home.uchicago.edu/~abney/abney_web/Software.html}
}
\seealso{
\code{\link[kinship2]{kinship}}, \code{\link[identity]{identity.coefs}}
}
\author{
Magnus Dehli Vigeland
}

\name{linkage.power}
\alias{linkage.power}
\alias{power.varyPar}
\title{
Power of a linkage study
}
\description{
Power analysis of parametric linkage studies
}
\usage{
linkage.power(x, N=100, all=FALSE, loop_breakers=NULL, threshold=NULL, seed=NULL)

%power.varyPar(x, N=100, varyPar, values, all=FALSE, loop_breakers=NULL, seed=NULL)
}

\arguments{
  \item{x}{
a \code{linkdat} object.
}
  \item{N}{
an integer; the number of markers to simulate.
}
  \item{all}{
a logical indicating whether markers should be simulated for all pedigree members. If FALSE, the availability vector \code{x$availability} will be used.
If \code{x} has existing marker data, the default is to simulate only the genotyped individuals.
}
  \item{loop_breakers}{
a numeric containing IDs of individuals to be used as loop breakers. Relevant only if the pedigree has loops. See \code{\link{breakLoops}}.
}
  \item{threshold}{
a numeric of length 1, or NULL. If not NULL, and if \code{maxonly=FALSE}, the output includes the percentage of simulated markers having LOD larger than \code{threshold}. 
}
  \item{seed}{
NULL, or a numeric seed for the random number generator. 
}
  % \item{varyPar}{
% a character indicating the parameter to vary. For autosomal models either 'f0', 'f1', 'f2', 'dfreq' or 'afreq'. For X-linked models
% either 'f0_m' (f0 for males), 'f1_m', 'f0_f', 'f1_f', 'f2_f', 'dfreq' or 'afreq'. 
% }
 % \item{values}{
% a numeric vector containing values for the varying parameter. 
% }
}
\value{
A \code{linkres} object, which is returned invisibly after a summary is printed

%The \code{power.varyPar} function creates a plot of the results and returns the maximum LOD score for each element of \code{values}.
}
\references{
The marker simulation takes place in \code{\link{SNPsim}}, which implements the algorithm used by SLINK: \url{http://linkage.rockefeller.edu/ott/SLINK.htm}
}
\author{
Magnus Dehli Vigeland
}

\seealso{
\code{\link{SNPsim}}, \code{\link{summary.linkres}}
}
\examples{
# Note: In the cases below, N is set very low. Increase to get more interesting results.
data(toyped)
x = linkdat(toyped, model=1)
linkage.power(x, N=5)

%power.varyPar(x, N=5, varyPar="f1", values=seq(0,1,length=11))

# X-linked example:
data(Xped)
y = linkdat(Xped, model=4)
linkage.power(y, N=1)

# Loop example:
z = addOffspring(cousinPed(1), father=7, mother=8, noffs=1, aff=2)
z = setModel(z, 2)
linkage.power(z, N=1, loop_breaker=7)
}
\keyword{ math }


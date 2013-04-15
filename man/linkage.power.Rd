\name{linkage.power}
\alias{linkage.power}
\alias{summary.powres}
\alias{power.varyPar}
\title{
Power of a linkage study
}
\description{
Power analysis of parametric linkage studies
}
\usage{
linkage.power(x, N=100, available=x$available, afreq=c(0.5, 0.5), 
              loop_breakers=NULL, threshold=NULL, seed=NULL, verbose=FALSE)

\method{summary}{powres}(object, threshold=NULL, \dots)
%power.varyPar(x, N=100, varyPar, values, all=FALSE, loop_breakers=NULL, seed=NULL)
}

\arguments{
  \item{x}{
a \code{\link{linkdat}} object with a valid \code{model}. (See \code{\link{setModel}}.) 
}
  \item{N}{
an integer; the number of markers to simulate.
}
  \item{available}{
a vector containing IDs of the available individuals, i.e. those whose genotypes should be simulated.
}
  \item{afreq}{
a numerical vector with sum 1; the population frequencies for the marker alleles.
}
  \item{loop_breakers}{
a numeric containing IDs of individuals to be used as loop breakers. Relevant only if the pedigree has loops. See \code{\link{breakLoops}}.
}
  \item{threshold}{
NULL, or a single numeric. If numeric, the output includes the percentage of simulated markers having LOD larger than \code{threshold}. 
}
  \item{seed}{
NULL, or a numeric seed for the random number generator. 
}
  \item{verbose}{
a logical passed on to \code{linkageSim}. If \code{TRUE}, some details are shown during the marker simulation. 
}
 \item{object}{
a \code{powres} object, normally produced by \code{linkage.power}. 
}
 \item{\dots}{
not used. 
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
The function prints a summary and returns invisibly a \code{powres} object, which is a list with the following entries:
\item{sim}{A \code{linkdat} object with the simulated markers}
\item{lod}{The LOD scores (computed with recombination fraction theta=0) of the simulated markers}
\item{maxlod}{The highest LOD score of the simulated markers}
\item{elod}{The average LOD score for the simulated markers} 

%The \code{power.varyPar} function creates a plot of the results and returns the maximum LOD score for each element of \code{values}.
}
\references{
The marker simulation takes place in \code{\link{linkageSim}}, which implements the algorithm used by SLINK: \url{http://linkage.rockefeller.edu/ott/SLINK.htm}
}
\author{
Magnus Dehli Vigeland
}

\seealso{
\code{\link{linkdat}}, \code{\link{linkageSim}}
}
\examples{
# Note: In the examples below N is set very low in order to reduce time consumption. Increase N to get more interesting results.
x = nuclearPed(3)
x = swapAff(x, c(1,3,4))
x = setModel(x, 1) # Autosomal dominant
linkage.power(x, N=1)

%power.varyPar(x, N=5, varyPar="f1", values=seq(0,1,length=11))

# X-linked recessive example:
data(Xped)
y = linkdat(Xped, model=4)
linkage.power(y, N=1)

# Power of homozygosity mapping:
z = addOffspring(cousinPed(1), father=7, mother=8, noffs=1, aff=2)
z = setModel(z, 2) # Autosomal recessive model
pow = linkage.power(z, N=1, loop_breaker=7, seed=123)
stopifnot(round(pow$maxlod, 1) == 1.2)
}
\keyword{ math }


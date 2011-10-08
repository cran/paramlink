\name{lod}
\alias{lod}
\title{
Two-point LOD score
}
\description{
Calculates the two-point LOD scores of a pedigree for the specified markers. The recombination ratio between the disease and marker loci can be either fixed at specific values, or optimized. 
}
\usage{
lod(x, markers=seq_len(x$nMark), t=c(0, 0.1, 0.2, 0.5), tol=0.01, 
    max.only=FALSE, silent=max.only)
}
\arguments{
  \item{x}{a \code{linkdat} object.
}
  \item{markers}{
an integer vector denoting which markers to use.
}
  \item{t}{
either a numeric containing specific recombination ratio(s), or the word "max", indicating that the recombination ratio should be optimized by the program.
}
	\item{tol}{
a numeric passed on to \code{\link{optimize}} as its tolerance parameter.
}
	\item{max.only}{
a logical indicating whether only the maximum LOD score should be returned.
}
	\item{silent}{
a logical. The output is less verbose if TRUE.
}
}
\details{
The LOD score is defined as LOD(t) = \code{log_10(L(t) - log_10(L(t=0.5)}, where \code{L(t)} denotes the likelihood of the observed pedigree data as a function of the recombination ratio \code{t} between the marker and the disease locus.
}
\value{
If \code{max.only=FALSE}, a \code{linkres} object, essentially a matrix containing the LOD scores. If \code{t} is numeric, the matrix has dimensions \code{length(t), length(markers)}, and the entry in row \code{t_0}, column \code{m_0} is the lod score of the pedigree for marker \code{m_0} assumming a recombination rate of \code{t_0}. 

If \code{t="max"}, the \code{linkres} matrix has one column per marker and two rows: The first containing the LOD score and the second the optimal recombination ratio for each marker.

If a marker has incompatible values (i.e. if a child of homozygous 1/1 parents has a 2 allele), the corresponding matrix entries are NaN.

If \code{max.only=TRUE}, only the highest LOD score is returned, as a numeric of length 1.
}
\author{
Magnus Dehli Vigeland
}
\seealso{
\code{\link{likelihood}}, \code{\link{optimize}}
}
\examples{
data(toyped)
x = linkdat(toyped, model=1)
res_max = lod(x, max.only=TRUE)
res_def = lod(x)[1,1]
res_opt = lod(x, t="max")[1,1]

stopifnot(all(0.3 == round(c(res_max, res_def, res_opt), 1)))

data(dominant)
y = linkdat(dominant)
y = setModel(y, model=1, penetrances=c(.001, .9, .99))
lod(y, markers=305:310)
lod(y, markers=305:310, t="max")
}
\keyword{ math }
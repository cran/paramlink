\name{lod}
\alias{lod}
\title{
Two-point LOD score
}
\description{
Calculates the two-point LOD scores of a pedigree for the specified markers. The recombination ratio between the disease and marker loci can be either fixed at specific values, or optimized. 
}
\usage{
lod(x, markers=seq_len(x$nMark), theta=0, loop_breakers=NULL, 
    max.only=FALSE, verbose=FALSE, tol=0.01)
}
\arguments{
  \item{x}{a \code{\link{linkdat}} object.
}
  \item{markers}{
an integer vector denoting which markers to use.
}
  \item{theta}{
either a numeric containing specific recombination ratio(s), or the word "max", indicating that the recombination ratio should be optimized by the program.
}
  \item{loop_breakers}{
a numeric containing IDs of individuals to be used as loop breakers. Relevant only if the pedigree has loops. See \code{\link{breakLoops}}.
}
  \item{max.only}{
a logical indicating whether only the maximum LOD score should be returned.
}
  \item{verbose}{
a logical: verbose output or not.
}
  \item{tol}{
a numeric passed on to \code{\link{optimize}} as its tolerance parameter.
}
}
\details{
The LOD score is defined as LOD(theta) = \code{log_10(L(theta) - log_10(L(theta=0.5)}, where \code{L(theta)} denotes the likelihood of the observed pedigree data as a function of the recombination ratio \code{theta} between the marker and the disease locus.
}
\value{
If \code{max.only=FALSE}, a \code{linkres} object, essentially a matrix containing the LOD scores. If \code{theta} is numeric, the matrix has dimensions \code{length(theta), length(markers)}, and the entry in row \code{t_0}, column \code{m_0} is the lod score of the pedigree for marker \code{m_0} assumming a recombination rate of \code{t_0}. 

If \code{theta="max"}, the \code{linkres} matrix has one column per marker and two rows: The first containing the LOD score and the second the optimal recombination ratio for each marker.

If a marker has incompatible values (i.e. if a child of homozygous 1/1 parents has a 2 allele), the corresponding matrix entries are NaN.

If \code{max.only=TRUE}, only the highest LOD score is returned, as a numeric of length 1.
}
\author{
Magnus Dehli Vigeland
}
\seealso{
\code{\link{likelihood}}, \code{\link{optimize}}, \code{\link{breakLoops}}
}
\examples{
data(toyped)
x = linkdat(toyped, model=1)
res = lod(x)
res_theta = lod(x, theta=c(0, 0.1, 0.2, 0.5))
res_max = lod(x, theta="max")
stopifnot(all(0.3 == round(c(res, res_theta['0',], res_max['LOD',]), 1)))

# bigger pedigree with several markers
data(dominant)
y = linkdat(dominant)
y = setModel(y, model=1, penetrances=c(.001, .9, .99))
lod(y, markers=305:310)
lod(y, markers=305:310, theta="max")

# Example with pedigree with loops:
data(twoloops)
z = linkdat(twoloops, model=2) # fully penetrant autosomal recessive model.
z = addMarker(z, 0, alleles=1:2, afreq=c(0.001, 0.999)) # add (empty) diallelic marker
z = modifyMarker(z, marker=1, ids=15:16, genotype=1) # make individuals 15 and 16 homozygous for the rare allele.
res1 = lod(z, loop_breakers=c(8,12))

z = modifyMarker(z, marker=1, ids=c(7,9,11,13), genotype=3, alleles=1:3, afreq=c(0.001, 0.499, 0.5)) # making the marker triallelic and adding some genotypes.
plot(z,1)
res2 = lod(z, loop_breakers=c(8,12))
z$markerdata[[1]]
z = modifyMarker(z, marker=1, alleles=1:4, afreq=c(0.001, 0.499, 0.25, 0.25)) 
z$markerdata[[1]]
res3 = lod(z, loop_breakers=c(8,12))
res1;res2;res3
stopifnot(all(3 == round(c(res1,res2,res3), 1)))
}
\keyword{ math }
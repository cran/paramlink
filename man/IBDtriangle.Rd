\name{IBDtriangle}
\alias{IBDtriangle}
\title{
IBD triangle plot
}
\description{
This function draws an IBD triangle, as typically used to visualize pairwise relatedness of non-inbred individuals. 
Various annotations are available, including points marking the most common relationships, contour lines for the kinship coefficients,
and shading of the unattainable region.
}
\usage{
IBDtriangle(relationships = c("UN", "PO", "MZ", "S", "H,U,G", "FC", "SC", "DFC", "Q"), 
            kinship.lines = numeric(), shading = "lightgray", pch = 16, 
            cex_points = 1.2, cex_text = 1, axes = FALSE)
}
\arguments{
  \item{relationships}{
A character vector indicating relationships points to be included in the plot. 
UN=unrelated; PO=parent/offspring; MZ=monozygotic twins; S=full siblings; H=half sibling; U=uncle/aunt;
}
  \item{kinship.lines}{
A numeric vector. 
}
  \item{shading}{
The shading color for the unattainable region.
}
  \item{pch}{
Symbol used for the relationship points. (See \code{\link{par}}.)
}
  \item{cex_points}{
A single numeric controlling the symbol size for the relationship points.
}
  \item{cex_text}{
A single numeric controlling the font size for the relationship labels.
}
  \item{axes}{
Draw surrounding axis box?
}
}
\details{
For any pair of non-inbred individuals A and B, their genetic relationship can be summarized by the IBD coefficients 
\eqn{(\kappa_0, \kappa_1, \kappa_2)}, where
\deqn{\kappa_i = P(A and B share i alleles IBD at random autosomal locus).}
Since \eqn{\kappa_0+\kappa_1+\kappa_2=1}, any relationship corresponds to a point in the triangle in the \eqn{(\kappa_0, \kappa_2)}-plane
defined by \eqn{\kappa_0 \ge 0, \kappa_2 \ge 0, \kappa_0 + \kappa_2 \le 1}. The choice of \eqn{\kappa_0} and \eqn{\kappa_2} as the axis variables
is done for reasons of symmetry and is not significant (other authors have used different views of the triangle).

As shown in Thompson (76) points in the subset of the triangle defined by \eqn{4\kappa_0\kappa_2 > \kappa_1} is unattainable for pairwise relationships. 
By default this region in shaded in a "lightgray" color.

The IBD coefficients are linearly related to the kinship coefficient \eqn{phi} by the formula 
\deqn{\phi = 0.25\kappa_1 + 0.5\kappa_2.}
By indicating values for \eqn{phi} in the \code{kinship.lines} argument, the corresponding contour lines are shown as dashed lines in the triangle plot.
}
\value{
NULL
}
\references{
Thompson, E. A. (1975): The estimation of pairwise relationships. Annals of Human Genetics 39.

Thompson, E. A. (1976): A restriction on the space of genetic relationships. Annals of Human Genetics 40.
}
\author{
Magnus Dehli Vigeland
}
\seealso{
\code{\link{examineKinships}}
}
\examples{
IBDtriangle()
}
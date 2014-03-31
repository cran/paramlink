\name{as.matrix.linkdat}
\alias{as.matrix.linkdat}
\alias{restore_linkdat}
\title{
linkdat to matrix conversion
}
\description{
\code{as.matrix} creates a matrix representing the given \code{linkdat} object in standard LINKAGE format, 
with marker annotations and other info attached as attributes.

\code{restore_linkdat} is the reverse of \code{as.matrix}.
}
\usage{
\method{as.matrix}{linkdat}(x, include.attrs = TRUE, ...)

restore_linkdat(x, attrs = NULL)
}
\arguments{
 \item{x}{
a \code{\link{linkdat}} object. In \code{restore_linkdat}: A numerical matrix in LINKAGE format.
}
 \item{include.attrs}{
a logical indicating if marker annotations and other info should be attached as attributes. See value.
}
 \item{attrs}{
a list containing marker annotations and other \code{linkdat} info compatible with \code{x}, in the format produced by \code{as.matrix}.
If NULL, the attributes of \code{x} itself are used.
}
 \item{\dots}{
not used.
}
} 
\details{
The way \code{\link{linkdat}} objects are created in paramlink, marker data are stored as a list of \code{marker} objects. Each of these is 
essentially a matrix with various attributes like allele frequencies, map info a.s.o.. This format works well for marker-by-marker
operations (e.g. likelihoods and LOD scores), but makes it somewhat awkward to operate "horizontally", i.e. individual-by-individual, 
for instance if one wants to delete all genotypes of a certain individual, or rearrange the pedigree in some way. 

It is therefore recommended to convert the \code{linkdat} object to a matrix first, do the necessary manipulations on the matrix, and 
finally use \code{restore_linkdat}. Attributes are often deleted during matrix manipulation, so it may be necessary to store them in a
variable and feed them manually to \code{restore_linkdat} using the \code{attrs} argument. 

With default parameters, \code{restore_linkdat(as.matrix(x))} should reproduce \code{x} exactly.
}

\value{
For \code{as.matrix}: A matrix with \code{x$nInd} rows and \code{6 + 2*x$nMark} columns. 
The 6 first columns describe the pedigree in LINKAGE format, 
and the remaining columns contain marker alleles, using the internal (numerical) allele coding and 0 for missing alleles.
If \code{include.attrs = TRUE} the matrix has the attributes \code{markerattr} (a list of marker annotations), 
\code{available} (the availability vector) and \code{model} (the disease model, if present). 

For \code{restore_linkdat}: A \code{linkdat} object.
}

\author{
Magnus Dehli Vigeland
}

\seealso{
\code{\link{linkdat}}, \code{\link{as.data.frame.linkdat}}
}
\examples{
x = linkdat(toyped, model=1)
y = restore_linkdat(as.matrix(x))
stopifnot(all.equal(x,y))

# If attributes are lost during matrix manipulation: Use the 'attrs' argument.
xmatr = as.matrix(x)
newmatr = xmatr[-4, ] # NB: attributes are lost here
z = restore_linkdat(newmatr, attrs = attributes(xmatr))

# should be the same as:
z2 = removeIndividuals(x, 4)
stopifnot(all.equal(z, z2))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ manip }
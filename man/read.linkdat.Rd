\name{read.linkdat}
\alias{read.linkdat}
\alias{write.linkdat}
\title{
Read/write 'linkdat' objects
}
\description{
Create 'linkdat' object from a pedigree file with simulation and/or marker information (e.g. pre-makeped format). 
}
\usage{
read.linkdat(file, header=FALSE, missing=0, model=NULL, ...)

write.linkdat(x, file="", col.names=FALSE)
}

\arguments{
	\item{file}{
the name of the file containing the pedigree. See Details for format information.
}
  \item{header}{
a logical value indicating whether the file contains column names. 
}
  \item{missing}{
a numeric or character of length 1 indicating the code for missing marker alleles.
}
  \item{model}{
optional: the linkage model to be passed onto \code{\link{setModel}} (see this for details).
}
  \item{x}{
a code{linkdat} object.
}
  \item{col.names}{
a logical value indicating if column names should be included in the file. 
}
  \item{\dots}{
further arguments to be passed to read.table.
}
}
\details{
The input to \code{read.linkdat} should be a file with numerical columns in the following order (standard LINKAGE format): 
1 Family id (optional),
2 Individual id,
3 Father id,
4 Mother id,
5 Sex,
6 Affection status,
7 First allele of first marker,
8 Second allele of first marker,
9 First allele of second marker,
a.s.o.

Only columns 2-6 are mandatory. The first column is interpreted as family id if its first two entries are equal.
}
\value{
For \code{read.linkdat}, a \code{\link{linkdat}} object.
}
\author{
Magnus Dehli Vigeland
}
\seealso{
\code{\link{linkdat}}
}
\examples{
data(toyped)
toy = linkdat(toyped)
write.linkdat(toy, "toy.ped", col.names=TRUE)
toy2 = read.linkdat("toy.ped", header=TRUE)
stopifnot(all.equal(toy, toy2))
}
\keyword{ IO }

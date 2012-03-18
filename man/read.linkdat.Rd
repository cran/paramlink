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
read.linkdat(file, header = FALSE, simcol = NA, missing=0, 
             model=NULL, ...)

write.linkdat(x, file="", col.names=FALSE, famid=TRUE)
}

\arguments{
	\item{file}{
the name of the file containing the pedigree. See Details for format information.
}
	\item{x}{
a code{linkdat} object.
}
  \item{header}{
a logical value indicating whether the file contains column names. If TRUE, the names of columns 1-5 (or 2-6 if famid=TRUE) must be ID, FID, MID, SEX and AFF in this order.
}
  \item{famid}{
a logical value indicating whether the first column of the output file should be a family id, consisting of all 1's. (This makes the 
file of classical LINKAGE format, hence compatible with other linkage software.) If \code{header=TRUE} the column name is FAMID.
}
  \item{simcol}{
an integer indicating the column (if present) containing simulation statuses. Not necessary if \code{header=TRUE} and the header of the 
simulation column is SIM.
}
  \item{missing}{
a numeric or character of length 1 indicating the code for missing marker alleles.
}
  \item{model}{
optional: the linkage model to be passed onto \code{\link{setModel}} (see this for details).
}
  \item{col.names}{
a logical value indicating if column names should be included in the file. 
}
  \item{\dots}{
further arguments to be passed to read.table.
}
}
\details{
The columns should be: Family id (optional); individual id; father id; mother id; sex; affection status; simulation status (optional); marker1_allel1; marker1_allel2; a.s.o.
If the file includes headers, these must be ('FAMID'), 'ID', 'FID', 'MID', 'SEX', 'AFF' and 'SIM'. Marker column names are ignored.

If the file has no headers, the program looks for a family id column by checking if the numbers in the first column are all equal.
If this is found to be true, a statement to this effect is printed and the column is removed.
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
toy=linkdat(toyped, model=1)
write.linkdat(toy, "toy.ped", col.names=TRUE)
toy2=read.linkdat("toy.ped", header=TRUE)
toy2
write.linkdat(toy, "toy_famid.ped") 
toy3=read.linkdat("toy_famid.ped")
toy3
}
\keyword{ IO }

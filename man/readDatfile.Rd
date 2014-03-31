\name{readDatfile}
\alias{readDatfile}
\title{
Read dat file in LINKAGE format
}
\description{
Converts dat files in LINKAGE format to dat/map/freq files in MERLIN format
}
\usage{
readDatfile(datfile, chrom, write_to = NULL)
}
\arguments{
  \item{datfile}{
character. The path to the dat file.
}
  \item{chrom}{
integer chromosome number (needed to create the MERLIN map).
}
  \item{write_to}{
a character prefix used for naming the output files, or NULL if no files should be written.
}
}
\value{
If \code{write_to} is NULL, a list of data.frames named \code{dat}, \code{map} and \code{freq}.
}
\author{
Magnus Dehli Vigeland
}
\examples{
# No example given.
}
\keyword{ IO }

\name{twoloops}
\alias{twoloops}
\docType{data}
\title{
A consanguineous pedigree
}
\description{
A consanguineous pedigree.
}
\usage{data(twoloops)}
\format{
A data frame describing the pedigree in standard LINKAGE format (i.e. the columns are: Family ID, individual ID, father, mother, sex, affection status).
}
\examples{
data(twoloops)
x = linkdat(twoloops)
plot(x)
}
\keyword{datasets}

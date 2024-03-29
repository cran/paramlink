% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/datasets.R
\docType{data}
\name{Xped}
\alias{Xped}
\title{Example pedigree with X-linked disease pattern.}
\format{
A data frame with 15 rows and 6 columns. See \code{\link{toyped}} for
  details about the format.
}
\usage{
Xped
}
\description{
A complex pedigree with an X-linked recessive disease pattern
}
\details{
The format is standard LINKAGE (pre-makeped) format.
}
\examples{

Xped

# Convert to a 'linkdat' object and set a recessive X-linked model:
x = linkdat(Xped, model=4)
summary(x)
}
\keyword{datasets}

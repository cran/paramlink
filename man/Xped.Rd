\name{Xped}
\alias{Xped}
\docType{data}
\title{
Example pedigree with X-linked disease pattern.
}
\description{
A complex pedigree consistent with an X-linked recessive disease.
}
\usage{Xped}
\format{
  A data frame with one row for each individual, and the following columns:
  \describe{
    \item{\code{FAMID}}{family id}
    \item{\code{ID}}{individual id}
    \item{\code{FID}}{father id (0 if founder)}
    \item{\code{MID}}{mother id (0 if founder)}
    \item{\code{SEX}}{sex: male=1, female=2}
    \item{\code{AFF}}{affection status: healthy=1, affected=2, unknown=0}
  }
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

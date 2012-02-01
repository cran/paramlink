\name{toyped}
\alias{toyped}
\docType{data}
\title{
Toy pedigree for linkage analysis
}
\description{
Toy pedigree with 4 individuals typed at 1 SNP marker. Individual 1 is missing one allele.
}
\usage{data(toyped)}
\format{
  A data frame with one row for each individual, and the following columns:
  \describe{
    \item{\code{FAMID}}{family id}
    \item{\code{ID}}{individual id}
    \item{\code{FID}}{father id (0 if founder)}
    \item{\code{MID}}{mother id (0 if founder)}
    \item{\code{SEX}}{sex: male=1, female=2}
    \item{\code{AFF}}{affection status: healthy=1, affected=2, unknown=0}
    \item{\code{M_A1}}{first marker, first allele}
    \item{\code{M_A2}}{first marker, second allele }
  }
}
\details{
The format is standard LINKAGE (pre-makeped) format.
}

\examples{
data(toyped)
toyped
}
\keyword{datasets}

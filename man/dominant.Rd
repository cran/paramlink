\name{dominant}
\alias{dominant}
\docType{data}
\title{
Example pedigree for linkage analysis: 23 individuals, 650 SNPs
}
\description{
Pedigree showing an autosomal dominant inheritance pattern. 23 individuals, of which 15 are genotyped with 650 SNP markers.
}
\usage{data(dominant)}
\format{
  A data frame with 23 observations on 1305 variables. The columns are ID, FID MID, SEX, AFF, M1_1, M1_2, ..., M650_1, M650_2. The first 5 columns describe the pedigree structure and affection status, while the final 1300 columns contain the marker alleles. All markers are SNPs, whose alleles are written as 1 and 2. Missing alleles are denoted by 0.
}

\examples{
data(dominant)
x = linkdat(dominant)
summary(x)
}
\keyword{datasets}

\name{largefam}
\alias{largefam}
\docType{data}
\title{
Example pedigree for linkage analysis: 23 individuals, 650 SNPs
}
\description{
Medical pedigree with autosomal dominant inheritance pattern. 23 individuals, 15 of which are genotyped at 650 SNP markers.
}
\usage{data(largefam)}
\format{
  A data frame with 23 observations on 1305 variables. The columns are ID, FID MID, SEX, AFF, M1_1, M1_2, ..., M650_1, M650_2. The first 5 columns describe the pedigree structure and affection status, while the final 1300 contain the marker information. All alleles are written numerically, with missing = 0. For example, for a SNP marker each allele is either 0 (if missing), 1 or 2.
 
}

\examples{
data(largefam)
}
\keyword{datasets}

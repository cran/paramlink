\name{exclusionPower}
\alias{exclusionPower}
\title{
Power of exclusion
}
\description{
Computes the power (of a single marker) of excluding a claimed relationship, given the true relationship. 
}
\usage{
exclusionPower(ped_claim, ped_true, ids, alleles, afreq = NULL, known_genotypes = list(), loop_breakers = NULL, Xchrom = FALSE)
}
\arguments{
  \item{ped_claim}{
a \code{\link{linkdat}} object, or a list of such (if disconnected), describing the claimed relationship.
}
  \item{ped_true}{
a \code{\link{linkdat}} object, or a list of such (if disconnected), describing the true relationship. Names must be consistent with \code{ped_claim}.
}
  \item{ids}{
individuals available for genotyping.
}
  \item{alleles}{
a numeric or character vector containing marker alleles names
}
  \item{afreq}{
a numerical vector with allele frequencies. An error is given if they don't sum to 1 (rounded to 3 decimals).
}
  \item{known_genotypes}{
list of triplets \code{(a, b, c)}, indicating that individual \code{a} has genotype \code{b/c}.
}
  \item{loop_breakers}{
a numeric containing IDs of individuals to be used as loop breakers. Relevant only if the pedigree has loops. See \code{\link{breakLoops}}.
}
  \item{Xchrom}{
a logical: Is the marker on the X chromosome?
}
}
\details{
This function computes the "Power of exclusion", as defined and discussed in "A general approach to power calculation for
relationship testing" (Egeland et al., submitted 2012).
}
\value{
A single numeric value.
}
\author{
Magnus Dehli Vigeland
}
\examples{
### Example from Egeland et al. (2012): 
### Two females claim to be mother and daughter. We compute the power (for various markers) 
### to reject this claim if they in reality are sisters.

mother_daughter = nuclearPed(1, sex = 2)
sisters = relabel(nuclearPed(2, sex = c(2, 2)), c(101, 102, 2, 3))

# Equifrequent SNP:
exclusionPower(ped_claim = mother_daughter, ped_true = sisters, 
               ids = c(2, 3), alleles = 2)
# SNP with MAF = 0.1:
exclusionPower(ped_claim = mother_daughter, ped_true = sisters, 
               ids = c(2, 3), alleles = 2, afreq=c(0.9, 0.1))
# Equifrequent tetra-allelic marker:
exclusionPower(ped_claim = mother_daughter, ped_true = sisters, 
               ids = c(2, 3), alleles = 4)
# Tetra-allelic marker with one major allele:
exclusionPower(ped_claim = mother_daughter, ped_true = sisters, 
               ids = c(2, 3), alleles = 4, afreq=c(0.7, 0.1, 0.1, 0.1))              

# How does the power change if the true pedigree is inbred?               
sisters_LOOP = addParents(sisters, 101, father = 201, mother = 202)
sisters_LOOP = addParents(sisters_LOOP, 102, father = 201, mother = 203)

# Equifrequent SNP:
exclusionPower(ped_claim = mother_daughter, ped_true = sisters_LOOP, loop = 101,
               ids = c(2, 3), alleles = 2)
# SNP with MAF = 0.1:
exclusionPower(ped_claim = mother_daughter, ped_true = sisters_LOOP, loop = 101,
               ids = c(2, 3), alleles = 2, afreq=c(0.9, 0.1))
# Equifrequent tetra-allelic marker:
exclusionPower(ped_claim = mother_daughter, ped_true = sisters_LOOP, loop = 101,
               ids = c(2, 3), alleles = 4)
# Tetra-allelic marker with one major allele:
exclusionPower(ped_claim = mother_daughter, ped_true = sisters_LOOP, loop = 101, 
               ids = c(2, 3), alleles = 4, afreq=c(0.7, 0.1, 0.1, 0.1))   
}
\keyword{math}

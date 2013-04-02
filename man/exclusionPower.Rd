\name{exclusionPower}
\alias{exclusionPower}
\title{
Power of exclusion
}
\description{
Computes the power (of a single marker) of excluding a claimed relationship, given the true relationship. 
}
\usage{
exclusionPower(ped_claim, ped_true, ids, alleles, afreq = NULL, known_genotypes = list(), loop_breakers = NULL, Xchrom = F, plot = T)
}
\arguments{
  \item{ped_claim}{
a \code{\link{linkdat}} object, or a list of several linkdat and/or singleton objects, describing the claimed relationship. If a list, the sets of ID labels must be disjoint, 
that is, all ID labels must be unique. 
}
  \item{ped_true}{
a \code{\link{linkdat}} object, or a list of several linkdat and/or singleton objects, describing the true relationship. ID labels must be consistent with \code{ped_claim}.
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
a numeric containing IDs of individuals to be used as loop breakers. Relevant only if any of the pedigrees has loops. See \code{\link{breakLoops}}.
}
  \item{Xchrom}{
a logical: Is the marker on the X chromosome?
}
  \item{plot}{
either a logical or the character "plot_only", controlling if a plot should be produced. If "plot_only", a plot is drawn, but no further computations 
are done (useful for reproducing the plot in computer-intensive applications)
}
}
\details{
This function computes the "Power of exclusion", as defined and discussed in "A general approach to power calculation for
relationship testing" (Egeland et al., submitted 2012).
}
\value{
A single numeric value. If \code{plot="plot_only"}, the function returns NULL after producing the plot.
}
\author{
Magnus Dehli Vigeland
}
\examples{
### A standard case paternity case:
### Compute the power of exclusion when the claimed father is in fact unrelated to the child.
claim = nuclearPed(noffs=1, sex=2) # Specifies individual 1 as the father of 3
true = list(singleton(id=1,sex=1), singleton(id=3, sex=2)) # Specifies 1 and 3 as unrelated
available = c(1, 3) # Individuals 1 and 3 are available for genotyping

# Equifrequent autosomal SNP:
PE1 = exclusionPower(claim, true, available, alleles = 2, afreq=c(0.5,0.5))
PE1
# If the child is known to have genotype 1/1:
PE2 = exclusionPower(claim, true, available, alleles = 2, afreq=c(0.5,0.5), known_genotypes=list(c(3,1,1)))
PE2
# Equifrequent SNP on the X chromosome:
PE3 = exclusionPower(claim, true, available, alleles = 2, afreq=c(0.5,0.5), Xchrom=TRUE)
PE3
stopifnot(PE1==0.125, PE2==0.25, PE3==0.25)
### Example from Egeland et al. (2012): 
### Two females claim to be mother and daughter. We compute the power (for various markers) 
### to reject this claim if they in reality are sisters.
mother_daughter = nuclearPed(1, sex = 2)
sisters = relabel(nuclearPed(2, sex = c(2, 2)), c(101, 102, 2, 3))

# Equifrequent SNP:
PE1 = exclusionPower(ped_claim = mother_daughter, ped_true = sisters, 
               ids = c(2, 3), alleles = 2)
PE1
# SNP with MAF = 0.1:
PE2 = exclusionPower(ped_claim = mother_daughter, ped_true = sisters, 
               ids = c(2, 3), alleles = 2, afreq=c(0.9, 0.1))
PE2
# Equifrequent tetra-allelic marker:
PE3 = exclusionPower(ped_claim = mother_daughter, ped_true = sisters, 
               ids = c(2, 3), alleles = 4)
PE3
# Tetra-allelic marker with one major allele:
PE4 = exclusionPower(ped_claim = mother_daughter, ped_true = sisters, 
               ids = c(2, 3), alleles = 4, afreq=c(0.7, 0.1, 0.1, 0.1))              

PE4
# How does the power change if the true pedigree is inbred?               
sisters_LOOP = addParents(sisters, 101, father = 201, mother = 202)
sisters_LOOP = addParents(sisters_LOOP, 102, father = 201, mother = 203)

# Equifrequent SNP:
PE5 = exclusionPower(ped_claim = mother_daughter, ped_true = sisters_LOOP, loop = 101,
               ids = c(2, 3), alleles = 2)
PE5
# SNP with MAF = 0.1:
PE6 = exclusionPower(ped_claim = mother_daughter, ped_true = sisters_LOOP, loop = 101,
               ids = c(2, 3), alleles = 2, afreq=c(0.9, 0.1))
PE6
# Equifrequent tetra-allelic marker:
PE7 = exclusionPower(ped_claim = mother_daughter, ped_true = sisters_LOOP, loop = 101,
               ids = c(2, 3), alleles = 4)
PE7
# Tetra-allelic marker with one major allele:
PE8 = exclusionPower(ped_claim = mother_daughter, ped_true = sisters_LOOP, loop = 101, 
               ids = c(2, 3), alleles = 4, afreq=c(0.7, 0.1, 0.1, 0.1))   
PE8
stopifnot(round(c(PE1,PE2,PE3,PE4,PE5,PE6,PE7,PE8), 5) == 
      c(0.03125, 0.00405, 0.08203, 0.03090, 0.03125, 0.00765, 0.07617, 0.03457))
}
\keyword{math}

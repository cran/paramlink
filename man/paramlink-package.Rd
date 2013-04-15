\name{paramlink-package}
\alias{paramlink-package}
\alias{paramlink}
\docType{package}
\title{
Parametric linkage analysis in R
}
\description{
Parametric linkage analysis and other likelihood-based pedigree analyses based on the Elston-Stewart algorithm, including singlepoint LOD scores, power analysis and computation of marker genotype distributions. Multipoint analysis is provided through a convenient MERLIN wrapper. Includes many utility functions for creating, manipulating and plotting pedigrees, markers and linkage results.
}
\details{
\tabular{ll}{
Package: \tab paramlink\cr
Type: \tab Package\cr
Version: \tab 0.9-1\cr
Date: \tab 2013-04-11\cr
Depends: \tab R (>= 3.0), kinship2\cr
SystemRequirements: \tab For the MERLIN wrapper, MERLIN (www.sph.umich.edu/csg/abecasis/merlin) must be installed and pointed to in the PATH environment variable.\cr
License: \tab GPL (>= 2)\cr
}

Index:
\preformatted{
Xped                    Example pedigree with X-linked disease pattern.
as.matrix.linkdat       linkdat to matrix conversion
as.data.frame.linkdat   linkdat to data.frame conversion
createPedigree          Create simple pedigrees
dominant                Example pedigree for linkage analysis: 23
                        individuals, 650 SNPs
exclusionPower          Power of exclusion
likelihood              Pedigree likelihood
linkage.power           Power of a linkage study
linkageSim              Simulate markers linked to a disease locus.
linkdat                 linkdat objects
linkdat.utils           linkdat utilities
linkres                 S3 methods for class 'linkres'.
lod                     Two-point LOD score
marker                  Marker functions
markerSim               Marker simulation
merlin                  MERLIN wrapper function
modifyPedigree          Modify the pedigree of 'linkdat' objects
oneMarkerDistribution   Genotype probability distribution
plot.linkdat            Plot pedigrees with genotypes
randomPed               Random pedigree
setAvailable            Functions for modifying availability vectors
setModel                Set, change or display the model parameters for
                        'linkdat' objects
simpleSim               Unconditional marker simulation
toyped                  Toy pedigree for linkage analysis
twoMarkerDistribution   Genotype probability distribution
twoloops                A consanguineous pedigree
}
}
\author{
Magnus Dehli Vigeland

Maintainer: Magnus Dehli Vigeland <magnusdv@studmed.uio.no>
}
\keyword{ package }
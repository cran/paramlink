\name{createPedigree}
\alias{createPedigree}
\alias{nuclearPed}
\alias{cousinPed}
\alias{halfCousinPed}
\title{
Create simple pedigrees
}
\description{
These are utility functions for creating rudimentary \code{linkdat} objects with common pedigree structures. 
The function \code{nuclearPed} creates a nuclear family with a specified number of male/female offspring, while \code{cousinPed} and \code{halfCousinPed} create pedigrees linking two (half) cousins of a specified degree.
}
\usage{
nuclearPed(noffs, sex)
cousinPed(degree)
halfCousinPed(degree)
}
\arguments{
  \item{noffs}{
a positive integer, the number of offspring in the nuclear family.
}
  \item{sex}{
a vector of length \code{noffs}; indicating the gender (1=male, 2=female) of the offspring. If missing, all offspring are taken to be males.
}
 \item{degree}{
a non-negative integer: 0=siblings, 1=first cousins; 2=second cousins, a.s.o.
}
}
\details{
All individuals are created as unaffected. Use \code{\link{swapAff}} to edit this (see Examples).
Use \code{\link{swapSex}} to change gender of pedigree members.
The cousins created by \code{cousinPed} and \code{halfCousinPed} have zero removal. Removals are easily contructed using \code{\link{addOffspring}} or \code{\link{removeIndividuals}} (see Examples). 
}
\note{
To make changes of sex/affection status, or to create cousin removals, it is often helpful to plot the pedigree.
}
\value{
a \code{\link{linkdat}} object.
}
\author{
Magnus Dehli Vigeland
}
\seealso{
\code{\link{swapAff}}, \code{\link{swapSex}}, \code{\link{removeIndividuals}}, \code{\link{addOffspring}}, \code{\link{relabel}}
}
\examples{
# A nuclear family with 2 boys and 3 girls, 
# where the father and the two boys are affected.
x = nuclearPed(noffs=5, sex=c(1,1,2,2,2))
x = swapAff(x, ids=c(1,3,4))

# A pedigree showing a pair of third cousins
y = cousinPed(3)
plot(y)

# Adding a child on one side, to create third cousins once removed. 
# The pedigree plot above helps identifying which individual 
# should have the extra child
y = addOffspring(y, father=15, noffs=1)

# Half sibs:
z = halfCousinPed(0)

# Alternative construction:
z2 = nuclearPed(noffs=1, sex=1)
z2 = addOffspring(z2, father=1, noffs=1, sex=2)
z2 = relabel(z2, new=1:5,old=c(1,2,4,3,5))
stopifnot(all.equal(z,z2))
}
\keyword{ manip }

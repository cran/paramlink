\name{createPedigree}
\alias{createPedigree}
\alias{nuclearPed}
\alias{cousinPed}
\alias{halfCousinPed}
\title{
Create simple pedigrees
}
\description{
These are utility functions for creating rudimentary 'linkdat' objects with common pedigree structures. 
The function \code{nuclearPed} creates a nuclear family with a specified number of male/female offspring, while \code{cousinPed} and \code{halfCousinPed} create pedigrees linking two (half) cousins of a specified degree.
}
\usage{
nuclearPed(boys, girls = 0)
cousinPed(degree)
halfCousinPed(degree)
}
\arguments{
  \item{boys}{
a nonnegative integer, the number of male offspring.
}
  \item{girls}{
a nonnegative integer, the number of female offspring.
}
 \item{degree}{
a non-negative integer: 0=siblings, 1=first cousins; 2=second cousins, a.s.o.
}
}
\details{
All individuals are created as unaffected. Use \code{\link{swapAff}} to edit this (see Examples).
Use \code{\link{swapSex}} to change gender of pedigree members.
The cousins created by \code{\link{cousinPed}} and \code{\link{halfCousinPed}} have zero removal. Removals are easily contructed using \code{\link{addOffspring}} or \code{\link{removeIndiv}} (see Examples). 
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
\code{\link{swapAff}}, \code{\link{swapSex}}, \code{\link{removeIndiv}}, \code{\link{addOffspring}}
}
\examples{
# A nuclear family with 2 boys and 3 girls, 
# where the father and the two boys are affected.
x = nuclearPed(boys=2, girls=3)
x2 = swapAff(x, ids=c(1,3,4))

# A pedigree showing a pair of third cousins
y = cousinPed(3)
plot(y)

# Adding a child on one side, to create third cousins once removed. 
# The pedigree plot above helps identifying which individual 
# should have the extra child
y2 = addOffspring(y, father=15, noffs=1)

# Half sibs:
z = halfCousinPed(0)
# Alternative construction:
z2 = nuclearPed(boys=1)
z2 = addOffspring(z2, father=1, noffs=1, sex=2)
}
\keyword{ manip }

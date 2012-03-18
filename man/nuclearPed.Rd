\name{nuclearPed}
\alias{nuclearPed}
\alias{cousinPed}
\title{
Create nuclear families or n-degree cousin pedigrees
}
\description{
These are utility functions for creating rudimentary 'linkdat' objects with common pedigree structures. 
The function \code{nuclearPed} creates a nuclear family with a specified number of boy/girl offspring, while cousinPed creates a pedigree
linking two cousins of specified degree.
}
\usage{
nuclearPed(boys, girls = 0)
cousinPed(degree=1)
}
\arguments{
  \item{boys}{
a nonnegative integer, the number of male offspring.
}
  \item{girls}{
a nonnegative integer, the number of female offspring.
}
 \item{degree}{
an positive integer: 1=first cousins; 2=second cousins, a.s.o.
}
}
\details{
All individuals are created as unaffected. Use \code{link{modifyPedigree}} to edit this (see Examples).

The function \code{cousinPed} creates pedigrees linking two male cousins with zero removal. 
Removals are easily contructed using \code{link{addChild}} or \code{link{removeIndiv}} (see Examples). 
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
\code{link{modifyPedigree}}, \code{link{removeIndiv}}, \code{link{plot.pedigree}}
}
\examples{
# A nuclear family with 2 boys and 3 girls, 
# where the father and the two boys are affected.
x=nuclearPed(boys=2, girls=3)
x=modifyPedigree(x, id=c(1,3,4), attrib='AFF')
plot(x)

# A pedigree showing a pair of third cousins
y=cousinPed(3)
plot(y)

# Adding a child on one side, to create third cousins once removed. 
# The pedigree plot above helps identifying which individual 
# should have the extra child
y = addChildren(y, father=16, mother=0)
plot(y)
}
\keyword{ manip }

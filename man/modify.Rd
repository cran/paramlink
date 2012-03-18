\name{modifyPedigree}
\alias{modifyPedigree}
\alias{addChildren}
\alias{removeIndiv}
\title{Modify 'linkdat' objects}
\description{
Functions to modify the pedigree data of a 'linkdat' object.
}
\usage{
modifyPedigree(x, id, attrib, newval=NULL)

addChildren(x, father, mother, children=1, sex=1, aff=1) 

removeIndiv(x, id)

}
\arguments{
  \item{x}{
a \code{\link{linkdat}} object
}
  \item{id}{
an integer vector, indicating the individual(s) to be modified.
}
  \item{attrib}{
the attribute to be changed, either "AFF" or "SEX" (lower case letters are allowed).
}
  \item{newval}{
the new values of the parameter/attribute. In \code{modifyParameter}, a numeric of length 1. 
In \code{modifyPedigree} either a numeric of length \code{length(id)}, or NULL (which results in switching 1 <-> 2 the values of \code{attrib}).
}
  \item{father, mother}{
integers indicating the ID's of parents. One of these (but not both) can be 0, in which case a new founder individual is created. 
}
   \item{children}{
an integer indicating the number of offspring to be created.
}
   \item{sex, aff}{
integer vectors indicating the gender and affection statuses of the offspring to be created (recycled if not of length \code{children}).
}
}
\details{
When removing an individual, all descendents are also removed as well as all founders remaining without offspring.
}
\value{
The modified \code{linkdat} object.
}
\author{
Magnus Dehli Vigeland
}

\seealso{
\code{\link{linkdat}}
}
\examples{
data(toyped)
x=linkdat(toyped, model=1)
summary(x)

#switching sex of both offspring:
x=modifyPedigree(x, id=c(3,4), attrib='SEX')

#adding two children (non-affected males) to individual 4
x=addChildren(x, father=0, mother=4, children=2)
x

#and removing one of them again 
x=removeIndiv(x, 7)
x

#adding an affected half-sib
x=addChildren(x, father=1, mother=0, sex=2, aff=2)
x

#the affected half-sib has genotype AB:
x=modifyMarker(x, 8, c(1,2))
x
plot(x,1)

lod(x)
}
\keyword{ manip }

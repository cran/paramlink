\name{modifyPedigree}
\alias{modifyPedigree}
\alias{swapSex}
\alias{swapAff}
\alias{addOffspring}
\alias{addParents}
\alias{removeIndividuals}
\alias{trim}
\alias{relabel}
\title{Modify the pedigree of 'linkdat' objects}
\description{
Functions to modify the pedigree of a 'linkdat' object.
}
\usage{
swapSex(x, ids)

swapAff(x, ids, newval=NULL)

addOffspring(x, father, mother, noffs, ids, sex=1, aff=1) 

addParents(x, id, father, mother)

removeIndividuals(x, ids, verbose=TRUE)

trim(x, keep="available", return.ids=FALSE)

relabel(x, new, old)
}
\arguments{
  \item{x}{
a \code{\link{linkdat}} object
}
  \item{id, ids}{
the individual(s) to be modified. In \code{addOffspring} the (optional) \code{ids} argument is used to specify ID labels for the offspring to be created.
}
  \item{newval}{
a numeric, indicating affection status values for the \code{ids} individuals: 1=unaffected, 2=affected, 0=unknown. If NULL, the affection statuses are swapped 1 <-> 2, hence the main use of the \code{newval} argument is to assign 0's.
}
  \item{father, mother}{
integers indicating the IDs of parents. If missing, a new founder individual is created (whose ID will be 1+the largest ID already in the pedigree).
}
 \item{noffs}{
a single integer indicating the number of offspring to be created.
}
   \item{sex, aff}{
integer vectors indicating the gender and affection statuses of the offspring to be created (recycled if less than \code{noffs} elements).
}
  \item{verbose}{
A logical: Verbose output or not.
}
 \item{keep}{
A character, either "available" (trimming the pedigree for unavailable members) or "affected" (trimming for unaffected members).
}
 \item{return.ids}{
A logical. If FALSE, the trimmed pedigree is returned as a new \code{linkdat} object. If TRUE, a vector containing the IDs of "removable" individuals is returned
}
  \item{new}{
a numeric containing new labels to replace those in \code{old}.
}
  \item{old}{
a numeric containing ID labels to be replaced by those in \code{new}. If missing, \code{old} is set to \code{x$orig.ids}, i.e. all members in their original order.
}
}
\details{
When removing an individual, all descendants are also removed as well as founders remaining without offspring.
}
\value{
The modified \code{linkdat} object.
}
\author{
Magnus Dehli Vigeland
}

\seealso{
\code{\link{linkdat}}, \code{\link{nuclearPed}}
}
\examples{
data(toyped)
x = linkdat(toyped)

# To see the effect of each command below, use plot(x) in between.
x = addParents(x, id=2, father=5, mother=6)

x = swapSex(x, c(1,5))
x = swapSex(x, c(2,6))

x = addOffspring(x, mother=6, noffs=2, id=c(7,10))
x = removeIndividuals(x, 3)
x = swapAff(x, c(4,10))

stopifnot(setequal(x$orig.ids, c(1,2,4,5,6,7,10,11)))

# Trimming a pedigree
data(dominant)
x = linkdat(dominant)
x_affectedOnly = trim(x, keep="affected")

unavail = trim(x, keep="available", return.ids=TRUE)
nonaff = trim(x, keep="affected", return.ids=TRUE)
stopifnot(setequal(unavail, c(5, 19:23)), setequal(nonaff, c(6:7, 12:13, 19:23)))
}
\keyword{ manip }

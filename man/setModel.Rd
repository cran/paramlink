\name{setModel}
\alias{setModel}
\alias{print.linkdat.model}

\title{
Set, change or display the model parameters for 'linkdat' objects
}
\description{
Functions to set, change and display model parameters involved in parametric linkage analysis.
}
\usage{
setModel(x, model = NULL, chrom = NULL, penetrances = NULL, dfreq = NULL)

\method{print}{linkdat.model}(x, ...)
}

\arguments{
  \item{x}{
in \code{setModel}: a \code{\link{linkdat}} object. In \code{print.linkdat.model}: a \code{linkdat.model} object. 
}
  \item{model}{
NULL, or an object of class \code{linkdat.model}, namely a list with elements 
\code{chrom}, \code{penetrances} and \code{dfreq}. 
In the \code{setModel} function, the \code{model} argument can be one of the integers 1-4, with the following meanings:

1 = autosomal dominant; fully penetrant, dfreq=1e-5

2 = autosomal recessive; fully penetrant, dfreq=1e-5

3 = X-linked dominant; fully penetrant, dfreq=1e-5

4 = X-linked recessive; fully penetrant, dfreq=1e-5
}
  \item{chrom}{
a character, either 'AUTOSOMAL' or 'X'. Lower case versions are allowed and will be converted automatically.
}
  \item{penetrances}{
if \code{chrom=='AUTOSOMAL'}: a numeric of length 3 - \code{(f0, f1, f2)} - 
where \code{fi} is the probability of being affected given \code{i} disease alleles. 

If \code{chrom=='X'}: a list of two vectors, containing the penetrances for each sex: 
\code{penetrances = list(male=c(f0, f1), female=c(f0, f1, f2))}. 
}
  \item{dfreq}{
the population frequency of the disease allele.
}
  \item{...}{further parameters}
}

\value{
\code{setModel} returns a new \code{linkdat} object, whose \code{model} entry is a \code{linkdat.model} object: A list
containing the given \code{chrom}, \code{penetrances} and \code{dfreq}. 
}

\author{
Magnus Dehli Vigeland
}
\seealso{
\code{\link{linkdat}}
}
\examples{
data(toyped)
x = linkdat(toyped)
x1 = setModel(x, model=1)
summary(x1)

# The shortcut 'model=1' above is equivalent to 
x2 = setModel(x, chrom="AUTOSOMAL", penetrances=c(0,1,1), dfreq=1e-5)
stopifnot(all.equal(x1, x2))

# X-linked recessive model:
y1 = setModel(x, model=4, dfreq=0.01)
summary(y1)

# Long version of the above:
y2 = setModel(x, chrom="X", penetrances=list(male=c(0,1), female=c(0,0,1)),
              dfreq=0.01)
stopifnot(all.equal(y1, y2))

stopifnot(all.equal(y1, setModel(x, y1$model)))
}
\keyword{ math }

\name{linkdat}
\alias{linkdat}
\alias{print.linkdat}
\alias{summary.linkdat}
\alias{as.data.frame.linkdat}
\title{
linkdat objects
}
\description{
Functions to create and display 'linkdat' objects.
}
\usage{
linkdat(x, model=NULL, verbose=TRUE, missing=0)

\method{print}{linkdat}(x, ..., markers)

\method{summary}{linkdat}(object, ...)

\method{as.data.frame}{linkdat}(x, ..., famid=F, markers=seq_len(x$nMark), alleles=NULL, 
        missing=NULL, singleCol=FALSE)

%\method{subset}{linkdat}(x, subset=1:x$nInd, ..., markers=seq_len(x$nMark))
}
\arguments{
 \item{x, object}{ a \code{linkdat} object, or a matrix/data.frame. In the latter case, \code{x} should contain containing columns corresponding to those indicated above. The order of the columns must be correct, as column names are ignored.
}
  \item{model}{
a \code{linkdat.model} object (typically \code{y$model} for some linkdat object \code{y}), or a single integer with the following meaning: 
1 = autosomal dominant; 2 = autosomal recessive; 3 = X-linked dominant; 4 = X-linked recessive.
In each of these cases, the disease allele frequency is set to 0.00001, and all markers are assumed to be diallelic and equifrequent.
}
 \item{verbose}{
a logical; verbose output or not.
}
  \item{markers}{
a numeric indicating which markers should be included.
}
  \item{famid}{
a logical indicating if the family identifier should be included as the first column.
}
  \item{alleles}{
a character containing the allele names, e.g. \code{alleles=c("A","B")}.
}
  \item{missing}{
a character of length 1 indicating how missing alleles should be represented.
}
  \item{singleCol}{
a logical: Should the two alleles for each marker be pasted into one column or kept in separate columns?
}
%  \item{subset}{
%a numeric containing the individuals in the sub-pedigree to be extracted. 
%}
  \item{\dots}{
further arguments.
}
}
\details{
A \code{linkdat} object is a list whose entries describe a single medical pedigree. Internally the individuals are relabeled as 1,2,..., but this should rarely be of concern to the end user. Some pedigree checking is done, but it is recommended to plot the pedigree before doing any analysis.
}

\value{
A \code{linkdat} object, which is essentially a list with some or all of the following entries:
 \item{pedigree }{\code{data.frame} describing the pedigree in standard linkage format, i.e. one row for each individual and 5 columns: ID, FID, MID, SEX, AFF. }
 \item{orig.ids}{the original individual id labels}
 \item{nInd}{the number of individuals in the pedigree}
 \item{founders}{vector of the founder individuals (NB: Internal labeling used)}
 \item{nonfounders}{vector of the nonfounder individuals (NB: Internal labeling used)}
 \item{subnucs}{list containing all (maximal) nuclear families in the pedigree. Each nuclear family is given as a vector of the form c(pivot, father, mother, child1, ...), where the pivot is either the id of the individual linking the nuclear family to the rest of the pedigree, or 0 if there are none. (NB: Internal labeling used)}
 \item{markerdata}{\code{matrix} containing the markerdata}
 \item{nMark}{the number of markers}
 \item{sim}{a numeric vector indicating the simulation status for each individual, using the standard SLINK codes: 2 if genotypes should be simulated for the individual, 0 if not. (This peculiar set of possible values stems from the original SLINK program, where the values 1 and 3 were used as well.)}
 \item{model}{a \code{linkdat.model} object, essentially a list containing the model parameters. See \code{\link{setModel}} for details}
 \item{initial_probs}{only if markers are diallelic (i.e. if \code{x$model$nallel==2}) - a list of length nInd containing "initial" 
 haplotype probabilities for each individual, given the penetrance values and allele frequencies}
 }

\author{
Magnus Dehli Vigeland
}

\seealso{
\code{\link{setModel}}, \code{\link{read.linkdat}}, \code{\link{write.linkdat}}
}
\examples{
data(toyped)
toyped
x=linkdat(toyped, model=1)
x
summary(x)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ IO }
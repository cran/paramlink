\name{merlin}
\alias{merlin}
\alias{merlinUnlikely}
\title{
MERLIN wrappers
}
\description{
These functions require MERLIN to be installed and correctly pointed to in the PATH environment variable. 
The required input files are automatically created (in the working directory) and deleted afterwards. The MERLIN command is executed via a call to system().  
}
\usage{
merlin(x, markers=seq_len(x$nMark), model=TRUE, theta=NULL, 
       options="", verbose=FALSE, generate.files=TRUE, 
       cleanup=generate.files, logfile="")

merlinUnlikely(x, remove=FALSE, verbose=!remove)
}

\arguments{
  \item{x}{
a \code{\link{linkdat}} object
}
  \item{markers}{
an integer vector indicating which markers to use (default: all).
} 
 \item{model}{
a logical: If TRUE (and x$model is not NULL), the file "merlin.model" is created and "--model merlin.model" is included to the MERLIN command.
}
  \item{theta}{
a numeric with values between 0 and 0.5: The recombination value(s) for which the LOD score is computed. 
The values of \code{theta} are converted to centiMorgan positions using the Haldane map function and included in the MERLIN command using the \code{--position} parameter. Works only for single markers (i.e. \code{markers} must consist of a single integer).
}
  \item{options}{
a character with additional options to the MERLIN command. See details.
}
  \item{verbose}{
a logical: Show MERLIN output and other information, or not.
}
  \item{generate.files}{
a logical. If TRUE, the files "merlin.ped", "merlin.dat", "merlin.map", "merlin.freq" and (if \code{model=TRUE}) "merlin.model" are created in the working directory.
}
  \item{cleanup}{
a logical: Should the files be deleted automatically after the MERLIN run?
}
  \item{logfile}{
a character. If this is given, the MERLIN screen output will be written to a file with this name.
}
  \item{remove}{
a logical. If FALSE, the function returns the indices of markers found to unlikely. 
If TRUE, a new \code{linkdat} object is returned, where the unlikely markers have been deleted.
}
}
\details{
For these functions to work, MERLIN must be installed and the path to merlin.exe included in the PATH variable. The \code{merlin} function is first and foremost a wrapper to the parametric linkage functionality of MERLIN. By default the following MERLIN command is run (via a call to \code{\link{system}}) after creating appropriate files in the current working directory:

\code{merlin -p _merlin.ped -d _merlin.dat -m _merlin.map -f _merlin.freq --model _merlin.model --tabulate --markerNames --quiet}

The resulting multipoint LOD scores are extracted from the output and returned in R as a \code{\link{linkres}} object. 

Additional command parameters can be passed on using the \code{options} argument (this is simply pasted onto the MERLIN command, so dashes must be included). For example, to obtain singlepoint LOD scores instead of multipoint, set \code{options="--singlepoint"}. (The singlepoint scores should agree with the results of \code{lod(x)}, except in cases where some individuals have partial genotypes (see Examples).)

If \code{model=FALSE} the \code{--model merlin.model} part is removed from the MERLIN command above. This is necessary for some calculations, e.g. likelihoods (see Examples).

The \code{merlinUnlikely} function is a wrapper for MERLIN's "--error" command. The syntax is similar to that of \code{\link{mendelianCheck}}.
}
\value{
If \code{model=TRUE}, a \code{\link{linkres}} object. Otherwise a character containing the complete MERLIN output.

For \code{merlinUnlikely}, a numeric containing the indices of the unlikely, or (if \code{remove=TRUE}) a new \code{linkdat} object where the unlikely markers are removed.
}
\references{
\url{http://csg.sph.umich.edu/abecasis/Merlin/}
}
\author{
Magnus Dehli Vigeland
}
\examples{
\dontrun{
x = linkdat(toyped, model=1)
x

# MERLIN treats partial genotypes (i.e. one known and one unknown allele) as missing:
lod_merlin = merlin(x)
lod_partial = lod(x)
x = modifyMarker(x, marker=1, ids=1, genotype=0)
lod_missing = lod(x)
stopifnot(lod_merlin == round(lod_missing, 3))

# Likelihood computation by MERLIN:
merlin(x, model=F, options="--lik")
}
}
\keyword{ manip }
lod <-
function(x, markers=seq_len(x$nMark), theta=0, loop_breakers=NULL, max.only=FALSE, verbose=FALSE, tol=0.01) {
	stopifnot(inherits(x, 'linkdat'))
   if(inherits(x, 'singleton')) stop("This function is not applicable to singleton objects.")
	if (is.null(x$model)) 	stop("No model set.")
	if (x$nMark==0) 		stop("No marker data indicated.")
	if (x$hasLoops)		{	
		if(is.null(loop_breakers)) stop("The pedigree has loops. Please indicate loop breakers.")
		x = breakLoops(x, loop_breakers)
	}
	if (min(markers)<1 || max(markers)>x$nMark)	stop("Nonexistent marker indicated.")

	nalleles = unlist(lapply(x$markerdata[markers], attr, 'nalleles'))
	if(any(nalleles > 2)) {
		if(is.character(theta) && theta=="max") stop("The option theta='max' is not implemented for markers with more than 2 alleles.")
		return(.lodm(x, markers=markers, theta=theta, loop_breakers=loop_breakers, max.only=max.only, verbose=verbose))
	}

	ilink <- function(x, marker) {
		log_denom <- stopl <- likelihoodSNP(x, marker, logbase=10, TR.MATR=.TRhalf)
		if (log_denom==-Inf) return(c(NaN, NaN))

		startl <- likelihoodSNP(x, marker, logbase=10, TR.MATR=.TRzero)
		optimal <- optimize(likelihoodSNP, c(0, 0.5), x=x, marker=marker, afreq=NULL, logbase=10, tol=tol, maximum=TRUE)
		if (optimal$objective > max(startl,stopl)) {
			log_numer <- optimal$objective; theta_max <- optimal$maximum
		} else {
		log_numer <- max(startl, stopl); theta_max <- c(0,.5)[which.max(c(startl, stopl))]
		}
		c(log_numer - log_denom, theta_max)
	}
	
	.TRzero <- .TRmatr(0, x$model$chrom); .TRhalf <- .TRmatr(0.5, x$model$chrom)
	
	map = .getMap(x, na.action=1, verbose=F)[markers, , drop=F]

	if (is.numeric(theta)) {
		stopifnot(max(theta)<=0.5, min(theta)>=0)
		if (verbose) 
			cat("Computing singlepoint LOD scores at each marker\nfor the following recombination ", ifelse(length(theta)==1, "value:\n", "values:\n"),
				paste("  theta =",theta,collapse="\n"), "\n", sep="")
		
		markerdata_list = x$markerdata[markers]
		trm_list = lapply(theta, .TRmatr, chrom=x$model$chrom)
		denoms = unlist(lapply(markerdata_list, function(m) likelihoodSNP(x, m, theta=NULL, logbase=10, TR.MATR=.TRhalf) ))
		numers = vapply(markerdata_list, function(m)  unlist(lapply(trm_list, function(TR) likelihoodSNP(x, marker=m, theta=NULL, logbase=10, TR.MATR=TR))) , FUN.VALUE=numeric(length(theta)))
		res = numers - rep(denoms, each=length(theta))
		res = structure(res, dim=c(length(theta), length(markers)), dimnames = list(theta, map$MARKER), analysis="mlink", map=map, class="linkres")
	} 
	else if (identical(theta,"max")) {
		if (verbose) cat("Computing singlepoint LOD scores for each marker,\nmaximizing over all recombination values.\n\n") 
		res = unlist(lapply(markers, ilink, x=x))
		res = structure(res, dim=c(2, length(markers)), dimnames = list(c("LOD", "t_max"), map$MARKER), analysis="ilink", map=map, class="linkres")
	}
	
	if (verbose) summary(res)
	if (max.only) ifelse(all(is.na(res)), NA, max(res, na.rm=TRUE))
	else res
}

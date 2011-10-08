lod <-
function(x, markers=seq_len(x$nMark), t=c(0, .1, .2, .5), tol=0.01, max.only=FALSE, silent=max.only) {
	stopifnot(class(x)=="linkdat")
	if (is.null(x$model)) 	stop("No model set.")
	if (x$nMark==0) 		stop("No marker data indicated.")
	if (x$model$nallel>2 || any(unlist(lapply(attr(x$markerdata, "alleles")[markers], length)) > 2)) 
							stop("Sorry - only diallelic markers allowed.")
	if (min(markers)<1 || max(markers)>x$nMark)
							stop("Nonexistent marker indicated.")
	
	ilink <- function(x, marker) {
		log_denom <- stopl <- likelihood(x, marker, logbase=10, TR.MATR=.TRhalf)
		if (log_denom==-Inf) return(c(NaN, NaN))

		startl <- likelihood(x, marker, logbase=10, TR.MATR=.TRzero)
		optimal <- optimize(likelihood, c(0, 0.5), x=x, marker=marker, logbase=10, tol=tol, maximum=TRUE)
		if (optimal$objective > max(startl,stopl)) {
			log_numer <- optimal$objective; theta_max <- optimal$maximum
		} else {
		log_numer <- max(startl, stopl); theta_max <- c(0,.5)[which.max(c(startl, stopl))]
		}
		c(log_numer - log_denom, theta_max)
	}
	
	.TRzero <- .TRmatr(0, x$model$chrom); .TRhalf <- .TRmatr(0.5, x$model$chrom)

	if (is.numeric(t)) {
		stopifnot(max(t)<=0.5, min(t)>=0)
		if (!silent) 
			cat("Computing singlepoint LOD scores at each marker\nfor the following recombination ", ifelse(length(t)==1, "value:\n", "values:\n"),
				paste("  t =",t,collapse="\n"), "\n", sep="")
		trm_list = lapply(t, .TRmatr, chrom=x$model$chrom)
		denoms = unlist(lapply(markers, likelihood, x=x, t=NULL, logbase=10, TR.MATR=.TRhalf))
		numers = vapply(markers, function(m) unlist(lapply(trm_list, likelihood, x=x, m=m, t=NULL, logbase=10)), FUN.VALUE=numeric(length(t)))
		res = numers - rep(denoms, each=length(t))
		res = structure(res, dim=c(length(t), length(markers)), dimnames = list(t, paste("M", markers, sep="")), analysis="mlink", class="linkres")
	} 
	else if (identical(t,"max")) {
		if (!silent) cat("Computing singlepoint LOD scores for each marker,\nmaximizing over all recombination values.\n\n") 
		res = sapply(markers, ilink, x=x)
		res = structure(res, dimnames = list(c("LOD", "t_max"), paste("M", markers,sep="")), analysis="ilink", class = "linkres")
	}
	
	if (!silent) summary(res)
	if (max.only) if(all(is.na(res))) return(NA) else return(max(res, na.rm=TRUE))
	invisible(res)
}


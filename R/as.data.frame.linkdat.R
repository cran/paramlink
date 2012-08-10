as.data.frame.linkdat <- function(x, ..., famid=F, markers=seq_len(x$nMark), alleles=NULL, missing=NULL, singleCol=FALSE, sep="") {
	p = relabel(x$pedigree, x$orig.ids)
	if(famid) p = cbind(FAMID=x$famid, p)
	
	genotypes <- matrix(numeric(0),nrow=x$nInd)
	if (!is.null(markers) && length(markers)>0) {
		if (min(markers) < 1 || max(markers) > x$nMark) stop("Invalid marker number(s)")
		chrom = ifelse(is.null(x$model), "AUTOSOMAL", x$model$chrom)		
			
		genotypes = .prettyMarkers(x$markerdata[markers], alleles=alleles, sep=sep, missing=missing, 
					singleCol=singleCol, sex=p[, 'SEX'])
	}
	data.frame(p, genotypes, stringsAsFactors = FALSE, ...)
}



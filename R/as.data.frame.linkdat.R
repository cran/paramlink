as.data.frame.linkdat <- function(x, ..., famid=F, markers=seq_len(x$nMark), alleles=NULL, missing=NULL, singleCol=FALSE) {
	p = relabel(x$pedigree, x$orig.ids)
	if(famid) p = cbind(FAMID=x$famid, p)
	
	genotypes <- matrix(numeric(0),nrow=x$nInd)
	if (!is.null(markers) && length(markers)>0) {
		if (min(markers) < 1 || max(markers) > x$nMark) stop("Invalid marker number(s)")
		m = x$markerdata
		if (is.null(alleles)) alleles=attr(m,"alleles")[markers]
		if (is.null(missing)) missing=attr(m,"missing")
		chrom = ifelse(is.null(x$model), "AUTOSOMAL", x$model$chrom)
		
		genotypes = .prettyMarkers(m[, sort.int(c(2*markers-1,2*markers))], alleles=alleles, 
									  missing=missing, mNames=paste('M', markers, sep=""), singleCol=singleCol,
									  chrom=chrom, sex=p[, 'SEX'])
	}
	return(data.frame(p, genotypes, stringsAsFactors = FALSE, ...))
}
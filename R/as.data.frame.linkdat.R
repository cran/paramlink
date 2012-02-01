as.data.frame.linkdat <- function(x, ..., famid=F, markers=seq_len(x$nMark), alleles=NULL, missing=NULL, singleCol=FALSE, sep="") {
	p = relabel(x$pedigree, x$orig.ids)
	if(famid) p = cbind(FAMID=x$famid, p)
	
	genotypes <- matrix(numeric(0),nrow=x$nInd)
	if (!is.null(markers) && length(markers)>0) {
		if (min(markers) < 1 || max(markers) > x$nMark) stop("Invalid marker number(s)")
		chrom = ifelse(is.null(x$model), "AUTOSOMAL", x$model$chrom)		
			
		genotypes = .prettyMarkers(x$markerdata[markers], alleles=alleles, sep=sep, missing=missing, 
					singleCol=singleCol, chrom=chrom, sex=p[, 'SEX'])
	}
	data.frame(p, genotypes, stringsAsFactors = FALSE, ...)
}


.as.annotated.matrix = function(x) {
	p = cbind(FAMID=x$famid, relabel(x$pedigree, x$orig.ids))
	if(x$nMark>0) {
		p = cbind(p, do.call(cbind, x$markerdata))
		attr(p, "markerattr") = lapply(x$markerdata, attributes)
	}
	attr(p, "available") = x$available
	attr(p, "model") = x$model
	p
}

.restore.linkdat = function(x, attrs=NULL) { #x = annotated matrix
	if(is.null(attrs)) attrs = attributes(x)
	y = linkdat(x[, 1:6], model = attrs$model, verbose=FALSE)
	markers = x[, -(1:6)]
	nMark = ncol(markers)/2
	if(nMark==0) y = setMarkers(y, NULL)
	else {
		markerattr = attrs$markerattr
		markerdata_list = lapply(seq_len(nMark), function(k) {
			m = markers[, c(2*k-1,2*k)]
			attributes(m) = c(markerattr[[k]][-1], list(dim=dim(m)))
			m
		})
		class(markerdata_list) = "markerdata"
		y = setMarkers(y, markerdata_list)
	}
	
	y$available = intersect(y$orig.ids, attrs$available)
	#y$map = attrs$map
	y
}


.as.pedigree = function(x, deceased=numeric(0), aff2=NULL) {
	ped = relabel(x$pedigree, x$orig.ids) 
	AFF = ped[, 'AFF']; AFF[AFF == 0] = NA; AFF = AFF - 1 #kinship2:::plot.pedigree uses aff codes 1=aff, 0=unaff, NA=unknown
	#if (all(AFF==0)) AFF = AFF+1 # because of a bug(?) in kinship2
	if(!is.null(aff2)) {
		stopifnot(is.numeric(aff2), length(aff2)==x$nInd, all(aff2 %in% 0:2))
		aff2[aff2==0] = NA; aff2 = aff2 - 1
		AFF = cbind(aff1 = AFF, aff2=aff2)
	}
	
	status = ifelse(x$orig.ids %in% deceased, 1, 0)
	pedigree(id=ped[,'ID'], dadid=ped[,'FID'], momid=ped[,'MID'], sex=ped[,'SEX'], affected=AFF, status=status)
}

plot.linkdat <- function(x, marker=NULL, alleles=NULL, sep="/", missing="-", 
	id.labels=x$orig.ids, title=paste("Family", x$famid), available=FALSE, col = 1,
	deceased=numeric(0), starred=numeric(0), aff2=NULL, margins=c(4.1, 1, 4.1, 1), ...) {
	
	pedigree = .as.pedigree(x, deceased=deceased, aff2=aff2)
	
	if(available) cols = ifelse(x$orig.ids %in% x$available, 2, 1)
	else cols = col
	
	strid = rep(id.labels, length.out=x$nInd)
	starred = .internalID(x, starred)
	strid[starred] = paste(strid[starred], "*", sep="")
	if (!is.null(marker)) {
		m = x$markerdata[marker]
		chrom = ifelse(is.null(x$model), 'AUTOSOMAL', x$model$chrom)
		gg = .prettyMarkers(m, alleles=alleles, sep=sep, missing=missing, singleCol=TRUE, chrom=chrom, sex=x$pedigree[, 'SEX'])
		geno = gg[,1]
		for(i in seq_along(marker)[-1]) geno = paste(geno, gg[, i], sep="\n")
		if(is.null(strid) || !any(nzchar(strid))) strid = geno else strid = paste(strid, geno, sep="\n")
	}	
	#par(xpd=T)
	plot(pedigree, id=strid, col=cols, mar=margins, ...)
	if(!is.null(title)) title(title)
}

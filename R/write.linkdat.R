write.linkdat=function(x, prefix="", what=c("ped", "map", "dat", "freq", "model"), merlin=FALSE) {
	generated.files = character(0)
	
	if(any(c("map", "dat", "freq") %in% what)) .map = .getMap(x, na.action=1, verbose=F)
	
	if("ped" %in% what) {
		if(merlin) alleles=1:100 else alleles = NULL #assuming no marker has more than 100 alleles
		.ped = as.data.frame(x, famid=TRUE, missing=0, alleles=alleles)
		write.table(.ped, pedname <- paste(prefix, "ped", sep="."), col.names=F, row.names=F, quote=F)
		generated.files = c(generated.files, pedname)
	}
	
	if("model" %in% what && !is.null(x$model)) {
		dfreq = format(x$model$dfreq, scientific=F, decimal.mark=".")
		penets = paste(format(x$model$penetrances, scientific=F, decimal.mark="."), collapse=",")
		.model = c("my_disease", dfreq, penets, "my_model")
		write(.model, modelname <- paste(prefix, "model", sep="."), sep=" \t", ncolumns=4)
		generated.files = c(generated.files, modelname)
	}
	
	if("map" %in% what) {
		write.table(.map, mapname <- paste(prefix, "map", sep="."), col.names=F, row.names=F, quote=F)
		generated.files = c(generated.files, mapname)
	}	
	
	if("dat" %in% what) {
		.dat = cbind(code=c("A", rep("M", nrow(.map))), value=c("my_disease", .map$MARKER))
		write.table(.dat, datname <- paste(prefix, "dat", sep="."), col.names=F, row.names=F, quote=F)
		generated.files = c(generated.files, datname)
	}	
	
	if("freq" %in% what) {
		col1=col2=col3=NULL
		for(i in seq_along(.map$MARKER)) {
			m = x$markerdata[[i]]
			nal = attr(m, 'nalleles')
			allele_labels = if (merlin) seq_len(nal) else attr(m, 'alleles')
			col1 = c(col1, "M", rep("A", nal))
			col2 = c(col2, .map$MARKER[i], allele_labels)
			col3 = c(col3, "", format(attr(m, 'afreq'), scientifit=F, digits=6))
		}
		.freq = cbind(col1, col2, col3)
		write.table(.freq, freqname <- paste(prefix, "freq", sep="."), col.names=F, row.names=F, quote=F)
		generated.files = c(generated.files, freqname)
	}
	invisible(generated.files)
}

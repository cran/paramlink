print.linkdat <- function(x, ..., markers) {
	if (missing(markers)) marker.nos = seq_len(min(x$nMark, 10)) 
	else {
		if(length(markers) > 0 && max(markers) > x$nMark) stop("Nonexisting marker(s) indicated") 
		marker.nos = markers
	}
	print(as.data.frame(x, markers=marker.nos, missing="-", singleCol=TRUE), ...)
	if (missing(markers) && x$nMark > 10)
		cat("\nOnly first 10 markers are shown. Use option 'markers=' to print specified markers.\n")

	if (is.null(x$model)) cat("\nModel parameters:\nNo model parameters set.\n")
	else print(x$model)
}

print.linkdat.model <- function(x, ...) {
	cat("\nModel parameters:\n")
	model = x
	switch(model$chrom, 
		AUTOSOMAL = cat("Autosomal inheritance with penetrances: (f0, f1, f2) =", paste('(',paste(model$penetrances, collapse=", "),')',sep=""), "\n"),
		X = cat("X-linked inheritance with penetrances:\n\tMales: (f0, f1) =", paste('(', paste(model$penetrances$male, collapse=", "),')',sep=""), 
				"\n\tFemales: (f0, f1, f2) =", paste('(', paste(model$penetrances$female, collapse=", "),')',sep=""), "\n")
	)
	cat("Disease allele frequency:", model$dfreq,"\n")
	cat("Marker allele frequencies:", model$afreq,"\n")
}

summary.linkdat <- function(object, ...) {
	x <- object
	cat("Pedigree:\n")
	cat(x$nInd,"individuals\n")
	cat(length(x$founders),"founders,",length(x$nonfounders),"nonfounders\n")
	if((ant<-length(x$subnucs)) > 0) cat(ant,"nuclear", ifelse(ant==1, "subfamily","subfamilies"),"\n")
	aff = x$pedigree[, 'AFF']
	if(all(aff==1)) cat("No pedigree members affected by disease\n")
	else cat(sum(aff==2), "affected by disease,", sum(aff==1), "unaffected,", sum(aff==0), "with unknown affection status\n") 
	
	cat("\nMarker data:\n", x$nMark, " markers in total\n", sep="")
	if (x$nMark > 0) {
		miss = which(rowSums(m <- x$markerdata)==0)
		cat(length(miss), "individuals with no available genotypes")
		if(length(miss)>0) {
			cat(":", paste(x$orig.ids[miss], collapse=", "), "\n")
			cat(round(sum(m[-miss,]==0)/length(m[-miss,])*100, 2),"% missing alleles (excluding ungenotyped individuals)\n")
		} else 
			cat("\n", sum(m==0)/length(m)*100, "% missing alleles\n", sep="")
	}
	cat("\nMarker simulation:\n")
	if (is.null(x$sim)) cat("No simulation indicated\n")
	else {
		s=(x$sim==2)
		if (all(s)) cat("All individuals marked for simulation\n")
		else cat("Individuals marked for simulation:", paste(x$orig.ids[s], collapse=", "), 
			   "\nIndividuals excluded from simulation:", paste(x$orig.ids[!s], collapse=", "), "\n")
	}
	
	if (is.null(x$model)) cat("\nModel parameters:\nNo model parameters set\n")
	else print(x$model)
}

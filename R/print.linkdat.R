print.linkdat <- function(x, ..., markers=seq_len(x$nMark)) {
	if (length(markers) > 0)
		if (max(markers) > x$nMark) stop("Nonexisting marker(s) indicated.") 
		
	df = as.data.frame(x, markers=markers, missing="-", singleCol=TRUE)
	print(df, ...)
	if (x$nMark > 0 && length(markers)==0)  
		cat("Markers not shown. Use option 'markers=' to print specified markers.\n")
	if (is.null(x$model)) cat("\nModel parameters:\nNo model parameters set\n")
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
	cat(ant<-length(x$subnucs),"nuclear", ifelse(ant==1, "subfamily","subfamilies"),"\n")
	
	cat("\nMarker data:\n", x$nMark, " markers in total\n", sep="")
	if (x$nMark > 0) {
		miss = which(rowSums(m <- x$markerdata)==0)
		cat(length(miss), "individuals with no available genotypes")
		if(length(miss)>0) {
			cat(": ID =", miss,"\n")
			cat(sum(m[-miss,]==0)/length(m[-miss,])*100,"% missing genotypes (excluding ungenotyped individuals)\n")
		} else 
			cat("\n", sum(m==0)/length(m)*100, "% missing genotypes\n", sep="")
	}
	cat("\nMarker simulation:\n")
	if (is.null(x$sim)) cat("No simulation indicated\n")
	else {
		s=(x$sim==2)
		if (sum(s)==x$nInd) cat("All individuals marked for simulation\n")
		else cat("Individuals marked for simulation:", which(s),"\nIndividuals not marked for simulation:", which(!s), "\n")
	}
	
	if (is.null(x$model)) cat("\nModel parameters:\nNo model parameters set\n")
	else print(x$model)
}

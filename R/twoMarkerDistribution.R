

twoMarkerDistribution <- function(x, id, partialmarker1, partialmarker2, theta, loop_breakers=NULL, eliminate=99, verbose=TRUE) {
	starttime = proc.time()
	if (!inherits(m1 <- partialmarker1, "marker"))
		if(is.numeric(m1) && length(m1)==1 && m1 <= x$nMark) m1 = x$markerdata[[m1]] else stop("The 'partialmarker1' must be a 'marker' object.")
	if (!inherits(m2 <- partialmarker2, "marker"))
		if(is.numeric(m2) && length(m2)==1 && m2 <= x$nMark) m2 = x$markerdata[[m2]] else stop("The 'partialmarker2' must be a 'marker' object.")
	
	x = setMarkers(x, structure(list(m1 <- partialmarker1, m2 <- partialmarker2), class="markerdata"))
	
	afreq1 = attr(m1, 'afreq'); alleles1 = attr(m1, 'alleles'); nall1 = attr(m1, 'nalleles')
	afreq2 = attr(m2, 'afreq'); alleles2 = attr(m2, 'alleles'); nall2 = attr(m2, 'nalleles')
	allgenos1 = .allGenotypes(nall1); allgenos2 = .allGenotypes(nall2)
	geno.names = list(paste(alleles1[allgenos1[, 1]], alleles1[allgenos1[, 2]], sep="/"), paste(alleles2[allgenos2[, 1]], alleles2[allgenos2[, 2]], sep="/"))
		
	if(verbose) {
		cat("Conditioning on the following marker data:\n")
		print(data.frame(ID=x$orig.ids, 
						M1=as.character(.prettyMarkers(list(m1), missing="-", singleCol=TRUE)),
						M2=as.character(.prettyMarkers(list(m2), missing="-", singleCol=TRUE))), row.names=FALSE)
		cat("\nAllele frequencies, marker 1:\n"); print(structure(afreq1, names=alleles1))
		cat("\nAllele frequencies, marker 2:\n"); print(structure(afreq2, names=alleles2))
		cat("\nRecombination rate between marker loci:", theta,"\n")
	}
		
	grid.subset = .my.grid(c(.make.grid.subset(x, m1, id, "AUTOSOMAL", make.grid=F), .make.grid.subset(x, m2, id, "AUTOSOMAL", make.grid=F))) #Do this before loop breaking, since eliminate2 works better WITH the loops.
	
	if (x$hasLoops)	{		
		if(is.null(lb <- loop_breakers)) 	 stop("The pedigree has loops. Please indicate loop breakers.")
		if(verbose) cat(ifelse(length(lb)==1, "\nBreaking loop at individual", "\nBreaking loops at individuals"), .prettycat(lb, "and"), "\n")
		x = breakLoops(x, lb)
		m1 = x$markerdata[[1]]; m2 = x$markerdata[[2]]
	}
	int.id = .internalID(x, id)
	
	marginal = likelihood(x, locus1=m1, locus2=m2, theta=theta, eliminate=eliminate)
	if(marginal==0) stop("Partial marker data is impossible")
	
	probs = array(0, dim=sapply(geno.names, length), dimnames=geno.names)
	probs[grid.subset] = apply(grid.subset, 1, function(allg_rows) {
		m1[int.id, ] = allgenos1[allg_rows[1], ]
		m2[int.id, ] = allgenos2[allg_rows[2], ]
		likelihood(x, locus1=m1, locus2=m2, theta=theta, eliminate=eliminate)
	})
	
	res = probs/marginal
	if(verbose) {
		cat("\nJoint genotype distribution at the two markers for individual ", id,":\n", sep="")
		print(round(res,4))
		cat("\nTotal time used: ", (proc.time() - starttime)[["elapsed"]], " seconds.\n", sep="")
		return(invisible(res))
	}
	else res
}


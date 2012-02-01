

oneMarkerDistribution <- function(x, ids, partialmarker, theta=0.5, loop_breakers=NULL, verbose=TRUE) {
	if (is.null(x$model)) if (theta < 0.5) stop("Model must be set first. Use setModel().") else x = setModel(x,1,penetrances=c(.5,.5,.5))
	if ((chrom <- x$model$chrom) != "AUTOSOMAL") stop("Only autosomal models implemented at the moment.")
	nids = length(ids)
	
	if (is.numeric(partialmarker) && length(partialmarker)==1) {
		if (partialmarker < 1 || partialmarker > x$nMark) stop("Indicated marker does not exist.")
		else	partialmarker = x$markerdata[[partialmarker]]
	}
	if (class(partialmarker) != "marker") stop("The 'partialmarker' must be either a single integer or a 'marker' object.")
	x = setMarkers(x, m <- partialmarker)
	afreq = attr(m, 'afreq')
	alleles = attr(m, "alleles")
	
	if(verbose) {
		cat("Assuming that the marker is autosomal.\n")
		if(theta < 0.5)	cat("Recombination rate between marker and disease locus:", theta,".\n", sep="")
		else if(any(x$pedigree[,'AFF'] != 1)) cat("Assuming theta=0.5, i.e that the marker is independent of disease status.\n")
		
		cat("\nPartial marker data:\n")
		print(data.frame(ID=x$orig.ids, GENO=.prettyMarkers(list(m), missing="-", singleCol=TRUE, chrom=chrom)), row.names=FALSE)
		cat("\nAllele frequencies:\n")
		print(structure(afreq, names=alleles))
	}
		
	if (x$hasLoops)	{		
		if(is.null(lb <- loop_breakers)) 	 stop("The pedigree has loops. Please indicate loop breakers.")
		if(verbose) cat(ifelse(length(lb)==1, "\nBreaking loop at individual", "\nBreaking loops at individuals"), .prettycat(lb, "and"), "\n")
		x = breakLoops(x, lb)
		m = x$markerdata[[1]]
	}
	
	int.ids = .internalID(x, ids)
	heterogenos = .comb2(nall <- attr(m, 'nalleles'))
	allgenos = rbind(cbind(1:nall,1:nall), heterogenos)
	genos = lapply(int.ids, function(id) {
				geno = m[id, ]
				if (all(geno == 0)) return(allgenos) #if both alleles are missing
				else if(all(geno != 0)) return(rbind(geno)) #if genotype is already known
				else return(cbind(geno[geno!=0], 1:nall)) #if one allele is missing
				})
	geno.names = lapply(genos, function(matr) paste(alleles[matr[,1]], alleles[matr[,2]], sep=""))
	indexgrid = .my.grid(lapply(genos, function(matr) 1:nrow(matr))) 
	
	probs = array(0, dim=sapply(geno.names, length), dimnames=geno.names)
	for (r in 1:nrow(indexgrid)) {
		indices = indexgrid[r,]
		for(i in 1:nids) m[int.ids[i], ] = genos[[i]][indices[i], ]

		start.dat = .likelihood_startdata(x, m)
		dat = start.dat
		probs[rbind(indices)] = .likelihood_multi(x, theta, dat=start.dat)
	}
	
	if (sum(probs)==0)
		warning("\nAll genotype probabilities zero. Mendelian error?") 
	res = probs/sum(probs)
	if(verbose) {
		cat(ifelse(nids==1, "\nGenotype probability distribution for individual ",
				"\nJoint genotype probability distribution for individuals "), .prettycat(ids,"and"), ":\n", sep="")
		print(round(res,4))
		return(invisible(res))
	}
	else res
}

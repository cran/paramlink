#Function to compute genotype probabilities of a single individual (id), 
#conditional on existing genotypes of other individuals.
#argument 'partialmarker' must be either NULL (indicating no given marker data), a single integer indicating an existing marker,
#or a matrix with two columns giving alleles for each individual (0 = unkwown) 
#The 't' argument gives recombination fraction between marker and disease. 
#The default value t=0.5 implies no linkage, i.e., affection status has no effect on the probabilities.

genoDistr <- function(x, id, partialmarker=NULL, t=0.5) {
	if (is.null(x$model)) stop("No model set.")
	
	mis=0
	if (is.null(partialmarker)) partialmarker = matrix(0, nrow=x$nInd, ncol=2)
	else if (is.numeric(partialmarker) && length(partialmarker)==1) {
		k = partialmarker
		if (k < 1 || k > x$nMark) stop("Indicated marker does not exist.")
		else {partialmarker = getMarkers(x,k); mis=attr(x$markerdata, "missing")}
	}
	else partialmarker = as.matrix(partialmarker)
	if (!identical(all.equal(dim(partialmarker),c(x$nInd, 2)), TRUE))
		stop("Wrong dimensions of 'partialmarker'. Must be either NULL, a single integer or a ", x$nInd, "* 2 matrix.")
	
	y = setMarkers(x, partialmarker, missing=mis)
	m = y$markerdata
	chrom = y$model$chrom
	sex = y$pedigree[, 'SEX']
	cat("Computing genotype probabilities for individual", id, "\nconditional on the following existing genotypes:\n")
	print(data.frame(ID=1:y$nInd, GENO=.prettyMarkers(m, missing="-", singleCol=TRUE, chrom=chrom, sex=sex)))
	if (t==0.5) 
		cat("\nUsing default recombination t=0.5, i.e. marker segregates independently of disease.\n")
		
	#recall genotype codes: 00 -> 0, 11 -> 1, 22 -> 2, 12/21 -> 3, 01/10 -> 4, 02/20 -> 5.
	partialgeno = .diallel2geno(m)
	genolist <- list(list('0'=1:2, 'A'=1, 'B'=2,NULL,NULL,NULL), list('0'=1:3, '1'=1, '2'=2, '3'=3, '4'=c(1,3), '5'=2:3))
	geno.compat = genolist[[(chrom=="AUTOSOMAL" || sex[id]==2) + 1]][[partialgeno[id] + 1]] #example: If 'id' has genotype 0/1 (code=4), the possible genotypes are 1 (=1/1) and 3 (=1/2).
	
	alleles = attr(m, "alleles")
	if (chrom=="AUTOSOMAL" || sex[id]==2) {
		all.geno.codes = 1:3 
		all.geno.names = paste(alleles[c(1,2,1)], alleles[c(1,2,2)], sep="")  #c(11, 22, 12)
		ordr = c(1,3,2)  #for output: change order from (11, 22, 12) to (11, 12, 22) 
	} else {
		all.geno.codes = 1:2 
		all.geno.names = alleles
		ordr = 1:2
	}
	probs <- sapply(all.geno.codes, function(g) { 
		if (g %in% geno.compat) {
			partialgeno[id]<-g   
			return(likelihood(y, singleNum.geno=partialgeno, t=t)) 
		}
		else return(0)
	})
	if (sum(probs)==0)
		stop("\nIndividual ", id, ": All genotype probabilities zero. Mendelian error?") 
	res = probs/sum(probs)
	names(res) = all.geno.names

	cat("\nResults:\n")
	res[ordr]
}
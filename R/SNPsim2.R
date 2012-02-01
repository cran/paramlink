.SNPsim2 <- function(x, N=1, partialmarker=NULL, loop_breakers=NULL, unique=FALSE, seed=NULL) {
	if (is.null(x$model)) {
		if (all(x$pedigree[,'AFF']==1)) {x=setModel(x,1); cat("Unaffected pedigree; assuming autosomal model.\n\n")}
		else stop("No model set. Use setModel().")
	}

	if (x$model$nallel>2) stop("Sorry - only diallelic markers allowed.")
	if (any(!is.numeric(N), length(N)>1, N%%1 != 0)) stop("N must be a positive integer.")
	
	if (is.null(partialmarker)) partialmarker = matrix(0, nrow=nInd, ncol=2) else if(ncol(partialmarker)!=2 || nrow(partialmarker)!=nInd) stop("Wrong dimensions of marker matrix.")
	x = setMarkers(x, partialmarker, missing=0); m = x$markerdata
	if(length(attr(m, "alleles")[[1]])>2) stop("Sorry - only diallelic markers allowed.")
	cat("Simulating SNPs completely linked to disease locus, \nconditional on the following genotypes:\n")
	print(data.frame(ID=x$orig.ids, GENO=.prettyMarkers(m, missing="-", singleCol=TRUE, chrom=chr, sex=x$pedigree[, 'SEX'])))
	
	if (loops <- x$hasLoops)	{		
		if(is.null(loop_breakers)) 	 stop("The pedigree has loops. Please indicate loop breakers.")
		x = breakLoops(x, loop_breakers)
	}
	nInd = x$nInd; chr = x$model$chrom; afreq = x$model$afreq
	if (is.null(x$sim)) sim=rep(2,nInd) else sim=x$sim
	if (!is.null(seed)) set.seed(seed)
	
	.TRzero <- .TRmatr(0, chr); .TRhalf <- .TRmatr(0.5, chr)
	zgeno = .diallel2geno(m)
	#simulate only indivs with i) no given genotype, ii) sim-status=2 or descendants with sim-status=2.
	sim_indivs <- seq_len(nInd)[sapply(seq_len(nInd), function(i) zgeno[i]==0 && (sim[i]==2 || any(sim[descendants(x, i, original.id=F)]==2)))]


	genoprobs <- function(x, partialgeno, id, values) {  
			#values are 1:3 for autosomal models, and either 1:2 (males) or 1:3 (females) for X-linked models
			#outputs vector of length |values| with genotype probs for indiv id given pedigree og partial genotype information
			probs <- sapply(values, function(g) { partialgeno[id]<-g;   likelihood(x, afreq=afreq, singleNum.geno=partialgeno, TR.MATR=.TRzero) } )
			if (sum(probs)==0) { print(cbind(x$pedigree, partialgeno)); stop("\nIndividual ",id,": All genotype probabilities zero. Mendelian error?") }
			probs
	}
	
	switch(chr,	
	AUTOSOMAL = {
		# simulation order: founders first (speeds up the likelihoods). 
		sim_indivs <- sim_indivs[order( sim_indivs %in% x$nonfounders)]
		
		# pre-calculate probabilities for the first 'init' individuals (big time saver!)
		init <- which.min(3^seq_along(sim_indivs) + 3*N*(length(sim_indivs)-seq_along(sim_indivs)))  #optimal 'init' minimizes the number of likelihood() calls
		initg <- t(.my.grid( rep(list(1:3), init ) ))
		initp <- apply(initg, 2, function(g) { 
			zgeno[ sim_indivs[1:init] ] <- g;   
			likelihood(x, afreq=afreq, singleNum.geno=zgeno, TR.MATR=.TRzero) } )
		if (identical(sum(initp), 0)) stop("All genotype probabilities zero. Wrong model?")
		
		# initialize marker matrix and pre-fill the rows of the 'init' individuals (i.e. sim_indivs[1:init]) 
		markers <- matrix(rep.int(zgeno, N), ncol=N)
		markers[sim_indivs[1:init], ] <- initg[ , suppressWarnings(sample.int( ncol(initg), size=N, replace=TRUE, prob=initp )) ]

		# do the rest of the individuals (only those present in sim_indivs)
		for (i in sim_indivs[-c(1:init)])
			markers[i, ] <- apply(markers, 2, function(partgeno) sample.int(3, size=1, prob=genoprobs(x, partgeno, id=i, values=1:3)) )
	},
	X = {
		ped = x$pedigree
		males = sim_indivs[ped[,'SEX'][sim_indivs]==1]; females = sim_indivs[ped[,'SEX'][sim_indivs]==2]
		males = males[order(males %in% x$nonfounders)]; females = females[order(females %in% x$nonfounders)]  #quicker with this?
		n_males=length(males); n_females=length(females)
		
		#find optimal 'init' values for males/females
		calls = outer(0:n_males, 0:n_females, function(m, f) 2^m * 3^f + 2*(n_males-m)*N + 3*(n_females-f)*N) # = number of times likelihood() is called.
		calls.min = arrayInd(which.min(calls), dim(calls)) 
		init_m = calls.min[1]-1; init_f = calls.min[2]-1
		init_indivs = c(males[seq_len(init_m)], females[seq_len(init_f)])
		#pre-calculate probabilities for the first 'init' individuals
		initg <- t(.my.grid( list(1:2, 1:3)[rep(1:2, c(init_m, init_f))] ))
		initp <- apply(initg, 2, function(g) { zgeno[ init_indivs ] <- g;   likelihood(x, afreq=afreq, singleNum.geno=zgeno, TR.MATR=.TRzero) } )

		#initialize marker matrix and pre-fill the rows of the 'init_indivs' 
		markers <- matrix(0, nrow=nInd, ncol=N)
		markers[init_indivs, ] <- initg[ , suppressWarnings(sample.int( ncol(initg), size=N, replace=TRUE, prob=initp)) ]

		#do the rest of the males (only those present in sim_indivs)
		for (i in setdiff(males, males[seq_len(init_m)])) 
			markers[i, ] <- apply(markers, 2, function(partgeno) sample.int(2, size=1, prob=genoprobs(x, partgeno, id=i, values=1:2)) )
		
		for (i in setdiff(females, females[seq_len(init_f)])) 
			markers[i, ] <- apply(markers, 2, function(partgeno) sample.int(3, size=1, prob=genoprobs(x, partgeno, id=i, values=1:3)) ) }
	)
	markers[(zgeno == 0) & (sim!=2), ] <- 0

	if (unique) markers <- unique(markers, MARGIN=2)
	x = setMarkers(x, .geno2diallel(markers))
	if(loops) x = tieLoops(x)
	x
}


SNPsim <- function(x, N=1, unique=FALSE, seed=NULL) {
	if (is.null(x$model)) stop("No model set.")
	if (x$model$nallel>2)
		stop("Sorry - only diallelic markers allowed.")
	if (any(!is.numeric(N), length(N)>1, N%%1 != 0)) stop("N must be a positive integer.")
	
	nInd = x$nInd; chr = x$model$chrom
	if (is.null(x$sim)) sim=rep(2,nInd) else sim=x$sim
	if (!is.null(seed)) set.seed(seed)
	
	#simulate only indivs who have sim-status=2 or have descendants with sim-status=2.
	#simulation order: founders first (speeds up the likelihoods).
	sim_indivs <- seq_len(nInd)[sapply(seq_len(nInd), function(i) sim[i]==2 || any(sim[descendants(x, i, original.id=F)]==2))]
	zgeno <- rep.int(0, nInd)
	.TRzero <- .TRmatr(0, chr); .TRhalf <- .TRmatr(0.5, chr)

	genoprobs <- function(x, partialgeno=rep.int(0,x$nInd), id, values) {  
			#values are 1:3 for autosomal models, and either 1:2 (males) or 1:3 (females) for X-linked models
			#outputs vector of length |values| with genotype probs for indiv id given pedigree og partial genotype information
			probs <- sapply(values, function(g) { partialgeno[id]<-g;   likelihood(x, singleNum.geno=partialgeno, TR.MATR=.TRzero) } )
			if (sum(probs)==0) { print(cbind(x$pedigree,partialgeno)); stop("\nIndividual ",id,": All genotype probabilities zero. Mendelian error?") }
			probs
	}
	
	switch(chr,	
	AUTOSOMAL = {
		sim_order <- sim_indivs[order( sim_indivs %in% x$nonfounders)]  #founders first.
		
		init = which.min(3^seq_along(sim_indivs) + 3*N*(length(sim_indivs)-seq_along(sim_indivs))) #gives the least number of times likelihood() calls.

		#pre-calculate probabilities for the first 'init' individuals
		initg <- t(expand.grid( rep(list(1:3), init ) ))
		initp <- apply(initg, 2, function(g) { zgeno[ sim_order[1:init] ] <- g;   likelihood(x, singleNum.geno=zgeno, TR.MATR=.TRzero) } )
		if (identical(sum(initp), 0)) stop("All genotype probabilities zero. Wrong model?")
		
		#initialize marker matrix and pre-fill the rows of the 'init' individuals (i.e. sim_order[1:init]) 
		markers <- matrix(0, nrow=nInd, ncol=N)
		markers[sim_order[1:init], ] <- initg[ , suppressWarnings(sample.int( ncol(initg), size=N, replace=TRUE, prob=initp )) ]

		#do the rest of the individuals (only those present in sim_order)
		for (i in sim_order[-c(1:init)])
			markers[i, ] <- apply(markers, 2, function(partgeno) sample.int(3, size=1, prob=genoprobs(x, partgeno, id=i, values=1:3)) ) }, 
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
		initg <- t(expand.grid( list(1:2, 1:3)[rep(1:2, c(init_m, init_f))] ))
		initp <- apply(initg, 2, function(g) { zgeno[ init_indivs ] <- g;   likelihood(x, singleNum.geno=zgeno, TR.MATR=.TRzero) } )

		#initialize marker matrix and pre-fill the rows of the 'init_indivs' 
		markers <- matrix(0, nrow=nInd, ncol=N)
		markers[init_indivs, ] <- initg[ , suppressWarnings(sample.int( ncol(initg), size=N, replace=TRUE, prob=initp)) ]

		#do the rest of the males (only those present in sim_indivs)
		for (i in setdiff(males, males[seq_len(init_m)])) 
			markers[i, ] <- apply(markers, 2, function(partgeno) sample.int(2, size=1, prob=genoprobs(x, partgeno, id=i, values=1:2)) )
		
		for (i in setdiff(females, females[seq_len(init_f)])) 
			markers[i, ] <- apply(markers, 2, function(partgeno) sample.int(3, size=1, prob=genoprobs(x, partgeno, id=i, values=1:3)) ) }
	)
	markers[sim!=2, ] <- 0
	if (unique) markers <- unique(markers, MARGIN=2)
	setMarkers(x, .geno2diallel(markers))
}
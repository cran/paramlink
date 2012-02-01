likelihood <-
function(x, marker, theta=NULL, afreq=NULL, logbase=NULL, TR.MATR=NULL, singleNum.geno=NULL) {

	.peel <- function(dat, sub, sexes, chrom, TR.MATR) { #dat = list(probs, haps)
		probs=dat[['probs']]; haps=dat[['haps']]
		fa_haps = haps[[sub['fa']]]; mo_haps = haps[[sub['mo']]]; fa_len = length(fa_haps); mo_len = length(mo_haps)
		likel <- probs[[ sub['fa'] ]] %*% .Internal(t.default( probs[[ sub['mo'] ]] ))  #outer product of parent probs
		piv=sub['pivot'] 
		piv_ind = .Internal(match(piv, sub[-1], nomatch=0, NULL)) 
		nonpiv.offs <- sub[-c(1:3, piv_ind+1)] #offspring except possible pivot
	
		switch(chrom, AUTOSOMAL = {
			for (b in nonpiv.offs) {
				mm = .Internal(aperm(	TR.MATR[ fa_haps, mo_haps, haps[[ b ]], drop=F ], c(3,1,2), TRUE)) * probs[[ b ]]
				likel <-  '*'(likel, .Internal(colSums(mm, length(haps[[b]]), fa_len*mo_len, FALSE) ))
			}	
			if (piv_ind==0) return(sum(likel))	
			else if (piv_ind==1) res = .Internal(rowSums(likel, fa_len, mo_len, FALSE))
			else if (piv_ind==2) res = .Internal(colSums(likel, fa_len, mo_len, FALSE))
			else {
				T = TR.MATR[ fa_haps, mo_haps, haps[[piv]], drop=F ]
				arr = .Internal(as.vector(T, 'any')) * .Internal(as.vector(likel, 'any')); dim(arr) = dim(T)
				res = .Internal(colSums(arr, fa_len*mo_len, length(haps[[piv]]), FALSE)) #sum for each entry of haps[[piv]]
				res = res * probs[[piv]]
			}
			haps[[piv]] <- haps[[piv]][res!=0];		probs[[piv]] <- res[res!=0]; 	return(dat=list(probs=probs, haps=haps))
		},
		X = {
			offlik = rep.int(1, length(mo_haps))
			for (moff in nonpiv.offs[ sexes[nonpiv.offs]==1 ])	{ #male offspring (except possible pivot)	
				mm = .Internal(t.default(TR.MATR[[1]][ mo_haps, haps[[moff]], drop=F ])) * probs[[moff]]
				offlik <- '*'(offlik, .Internal(colSums(mm, length(haps[[moff]]), mo_len, FALSE)) )
			}
			likel <- .Internal(t.default( .Internal(t.default(likel)) * offlik))  #multiply each ROW by final offlik
		
			for (foff in nonpiv.offs[ sexes[nonpiv.offs]==2 ]) { #female offspring (except possible pivot)
				mm = .Internal(aperm(	TR.MATR[[2]][ fa_haps, mo_haps, haps[[ foff ]], drop=F ], c(3,1,2), TRUE)) * probs[[ foff ]]
				likel <-  '*'(likel, .Internal(colSums(mm, length(haps[[foff]]), fa_len*mo_len, FALSE) ))
			}
			if (piv_ind==0) return(sum(likel))	
			else if (piv_ind==1) res = .Internal(rowSums(likel, fa_len, mo_len, FALSE))
			else if (piv_ind==2) res = .Internal(colSums(likel, fa_len, mo_len, FALSE))
			else 	#if pivot is offspring:
				if (sexes[piv]==1) { 		#if pivot is a male offspring
					T = TR.MATR[[1]][ mo_haps, haps[[piv]], drop=F ]
					res = .Internal(colSums(T * .Internal(colSums(likel, fa_len, mo_len, FALSE)), mo_len, length(haps[[piv]]), FALSE)) * probs[[piv]]   #res = colSums(T * colSums(likel))
				} else { 					#if pivot is a female offspring
					T = TR.MATR[[2]][ fa_haps, mo_haps, haps[[piv]], drop=F ]
					arr = .Internal(as.vector(T, 'any')) * .Internal(as.vector(likel, 'any')); dim(arr) = dim(T)
					res = .Internal(colSums(arr, fa_len*mo_len, length(haps[[piv]]), FALSE)) #sum for each entry of haps[[piv]]
					res = res * probs[[piv]]
				} 
			haps[[piv]] <- haps[[piv]][res!=0];		probs[[piv]] <- res[res!=0]; 	return(dat=list(probs=probs, haps=haps))
		})
	}
	
	if(x$hasLoops) stop("Unbroken loops in pedigree.")
	nInd=x$nInd; ped=x$pedigree; chr=x$model$chrom
	
	if (is.null(singleNum.geno)) {
      if (length(marker) == 1) marker = x$markerdata[[marker]]
		if (max(marker) > 2) stop("Marker has more than two alleles. You should use 'likelihood_multiallele' instead.")
		if (is.null(afreq)) afreq = attr(marker, 'afreq')
		singleNum.geno = .diallel2geno(marker)  #Coding genotypes as single integer: 00 -> 0, 11 -> 1, 22 -> 2, 12/21 -> 3, 01/10 -> 4, 02/20 -> 5.
	}
	if (is.null(afreq)) stop("Allele frequencies missing.")
	if (length(afreq)==1) afreq = c(afreq, 1- afreq)
	if (is.null(TR.MATR)) TR.MATR=.TRmatr(theta, chr)
	
	switch(chr, AUTOSOMAL = {		
		hap_list <- list('??'=1:10, 'AA'=1:3, 'BB'=8:10, 'AB'=4:7, 'A?'=1:7, 'B?'=4:10)[singleNum.geno + 1]  #+1 because coding starts at 0 (which is convenient in SNPsim a.s.o.)
		a=afreq[1]; b=afreq[2];
		init_p = x$initial_probs; 
		init_p[, x$founders] = init_p[, x$founders] * rep(c(a^2, 2*a*b, b^2), c(3, 4, 3))
		prob_list <- lapply(seq_len(nInd), function(i) init_p[, i][hap_list[[i]]])
	}, X = {	
		haplo.poss_X <- list(list('?'=1:4, 'A'=1:2, 'B'=3:4, NULL, NULL, NULL), list('??'=1:10, 'AA'=1:3, 'BB'=8:10, 'AB'=4:7, 'A?'=1:7, 'B?'=4:10))
		hap_list <- lapply(seq_len(nInd), function(i) haplo.poss_X[[ ped[i,'SEX'] ]][[singleNum.geno[i] + 1]])
		a=afreq[1]; b=afreq[2] 
		AfreqX <- list(male=c(a, a, b, b), female=rep(c(a^2, 2*a*b, b^2), c(3, 4, 3)))
		init_p_list = x$initial_probs
		for (i in x$founders) 	init_p_list[[i]] <- init_p_list[[i]] * AfreqX[[ ped[i,'SEX'] ]]
		prob_list <- lapply(seq_len(nInd), function(i) 	init_p_list[[i]][ hap_list[[i]] ])  #could use shorter mapply construction, but this is faster
	})
	hap_list <- lapply(seq_len(nInd), function(i) hap_list[[i]][prob_list[[i]]!=0])
	prob_list <- lapply(prob_list, function(v) v[v!=0])

	if(is.null(dups <- x$loop_breakers)) {
		dat = list(probs=prob_list, haps=hap_list)
		for (sub in x[['subnucs']]) 		
			dat = .peel(dat, sub, sexes=ped[, 'SEX'], chrom=chr, TR.MATR)
		likelihood = dat
	}
	else {
		origs = match(dups[, 1], x$orig.ids); copies = match(dups[, 2], x$orig.ids)
		sumover = .my.grid(lapply(seq_along(origs), function(i) intersect(hap_list[[origs[i]]], hap_list[[copies[i]]])))
		sumoverlist = lapply(seq_len(nrow(sumover)), function(ri) sumover[ri,])
		likelihood = 0
		for(r in sumoverlist) { 
			probs = prob_list; haps = hap_list
			for(i in seq_along(origs)) {
				haps[[origs[i]]] <- haps[[copies[i]]] <- r[i]
				probs[[origs[i]]] = probs[[origs[i]]][ r[i] == hap_list[[origs[i]]] ] #RIKTIG?? hap_list i stedet?
				probs[copies[i]] = list(1)
			}
			dat = list(probs=probs, haps=haps)
			for (sub in x[['subnucs']]) 	
				dat = .peel(dat, sub, sexes=ped[, 'SEX'], chrom=chr, TR.MATR)
			likelihood = likelihood + dat
		}
	}
	if (is.numeric(logbase)) log(likelihood, logbase) else likelihood
}

.TRmatr=function(theta, chrom) {
	if (is.null(theta)) stop("Argument 'theta' cannot be NULL.")
	haplo.single <- c('AD','AN','BD','BN')
	haplo.allpairs <- c('AADD','AADN','AANN','ABDD','ABDN','ABND','ABNN','BBDD','BBDN','BBNN')
	h <- c( c(1,.5,0,.5,.5*(1-theta),.5*theta,0,0,0,0), c(0,.5,1,0,.5*theta,.5*(1-theta),.5,0,0,0), c(0,0,0,.5,.5*theta,.5*(1-theta),0,1,.5,0), 
			c(0,0,0,0,.5*(1-theta),.5*theta,.5,0,.5,1) )
	dim(h) <- c(10,4); dimnames(h) <- list(haplo.allpairs, haplo.single)
	switch(chrom, 
	AUTOSOMAL = {
		T <- numeric(1000); dim(T) <- c(10,10,10); dimnames(T) <- list(haplo.allpairs, haplo.allpairs, haplo.allpairs)
		T[,,'AADD'] <- h[,'AD'] %*% t(h[,'AD'])
		T[,,'AADN'] <- h[,'AD'] %*% t(h[,'AN']) + h[,'AN'] %*% t(h[,'AD'])
		T[,,'AANN'] <- h[,'AN'] %*% t(h[,'AN'])
		T[,,'ABDD'] <- h[,'AD'] %*% t(h[,'BD']) + h[,'BD'] %*% t(h[,'AD'])
		T[,,'ABDN'] <- h[,'AD'] %*% t(h[,'BN']) + h[,'BN'] %*% t(h[,'AD'])
		T[,,'ABND'] <- h[,'AN'] %*% t(h[,'BD']) + h[,'BD'] %*% t(h[,'AN'])
		T[,,'ABNN'] <- h[,'AN'] %*% t(h[,'BN']) + h[,'BN'] %*% t(h[,'AN'])
		T[,,'BBDD'] <- h[,'BD'] %*% t(h[,'BD'])
		T[,,'BBDN'] <- h[,'BD'] %*% t(h[,'BN']) + h[,'BN'] %*% t(h[,'BD'])
		T[,,'BBNN'] <- h[,'BN'] %*% t(h[,'BN'])
		return(T)
	}, X = {
		TR_f <- numeric(400); dim(TR_f) <- c(4,10,10); dimnames(TR_f) <- list(haplo.single, haplo.allpairs, haplo.allpairs)
		TR_f['AD', , c('AADD','AADN','ABDD','ABDN')] <- h
		TR_f['AN', , c('AADN','AANN','ABND','ABNN')] <- h 
		TR_f['BD', , c('ABDD','ABND','BBDD','BBDN')] <- h 
		TR_f['BN', , c('ABDN','ABNN','BBDN','BBNN')] <- h
		return(list(h, TR_f))
	})
}

.my.grid = function (argslist, as.list=FALSE) {
	nargs <- length(argslist) 
    if (nargs == 0L) return(matrix(ncol=0, nrow=0))
    
    rep.fac <- 1L
    orep <- nr <- prod(unlist(lapply(argslist, length)))
	if(nr==0) return(matrix(ncol=nargs, nrow=0))
	
	res <- NULL
	for (x in argslist) {
		nx <- length(x)
		orep <- orep/nx
		res <- c(res, x[rep.int(rep.int(seq_len(nx), rep.int(rep.fac, nx)), orep)]) #this is res[, i]
		rep.fac <- rep.fac * nx
	}
	dim(res) <- c(nr, nargs)
	if(as.list) res = lapply(seq_len(nr), function(r) res[r,])
	res
}

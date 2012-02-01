

twoMarkerDistribution <- function(x, id, partialmarker1, partialmarker2, theta, loop_breakers=NULL, verbose=TRUE) {
	
	if (is.numeric(partialmarker1) && length(partialmarker1)==1) {
		if (partialmarker1 < 1 || partialmarker1 > x$nMark) stop("Indicated marker does not exist.")
		else	partialmarker1 = x$markerdata[[partialmarker1]]
	}
	
	if (is.numeric(partialmarker2) && length(partialmarker2)==1) {
		if (partialmarker2 < 1 || partialmarker2 > x$nMark) stop("Indicated marker does not exist.")
		else	partialmarker2 = x$markerdata[[partialmarker2]]
	}
	if (class(partialmarker1) != "marker") stop("'partialmarker1' must be either a single integer or a 'marker' object.")
	if (class(partialmarker2) != "marker") stop("'partialmarker2' must be either a single integer or a 'marker' object.")
	x = setMarkers(x, structure(list(m1 <- partialmarker1, m2 <- partialmarker2), class="markerdata"))
	
	afreq1 = attr(m1, 'afreq')
	afreq2 = attr(m2, 'afreq')

	if(verbose) {
		cat("Conditioning on the following marker data:\n")
		print(data.frame(ID=x$orig.ids, 
						M1=as.character(.prettyMarkers(list(m1), missing="-", singleCol=TRUE, chrom="AUTOSOMAL")),
						M2=as.character(.prettyMarkers(list(m2), missing="-", singleCol=TRUE, chrom="AUTOSOMAL"))), row.names=FALSE)
		cat("\nAllele frequencies, marker 1:\n"); print(structure(afreq1, names=attr(m1, 'alleles')))
		cat("\nAllele frequencies, marker 2:\n"); print(structure(afreq2, names=attr(m2, 'alleles')))
		cat("\nRecombination rate between marker loci:", theta,"\n")
	}
		
	if (x$hasLoops)	{		
		if(is.null(lb <- loop_breakers)) 	 stop("The pedigree has loops. Please indicate loop breakers.")
		if(verbose) cat(ifelse(length(lb)==1, "\nBreaking loop at individual", "\nBreaking loops at individuals"), .prettycat(lb, "and"), "\n")
		x = breakLoops(x, lb)
		m1 = x$markerdata[[1]]; m2 = x$markerdata[[2]]
	}
	
	int.id = .internalID(x, id)
	
	genos = lapply(list(m1,m2), function(mi) {
		alleles = attr(mi, 'alleles')
		nall = attr(mi, 'nalleles')
		geno = mi[int.id, ]
		res = switch(sum(geno != 0) + 1,  
			rbind(cbind(1:nall,1:nall), .comb2(nall)), #if both alleles are missing
			cbind(geno[geno!=0], 1:nall), #if one allele is missing
			rbind(geno) #if genotype is already known		
		)
		attr(res, 'gnames') = paste(alleles[res[,1]], alleles[res[,2]], sep="")
		res
	})
	
	len1 = nrow(genos[[1]]); len2 = nrow(genos[[2]])
	probs = matrix(0, nrow=len1, ncol=len2, dimnames=lapply(genos, attr, 'gnames'))

	for(g1 in 1:len1) for(g2 in 1:len2) {
		m1[int.id, ] = genos[[1]][g1,]
		m2[int.id, ] = genos[[2]][g2,]
		start.dat = .twopoint_general_startdata(x, m1, m2)
		probs[g1, g2] = .likelihood_twopoint_general(x, theta, dat=start.dat)	
	}
	
	if (sum(probs)==0)
		stop("\nIndividual ", id, ": All probabilities zero. Mendelian error?") 
	res = probs/sum(probs)
	if(verbose) {
		cat("\nJoint genotype distribution at the two markers for individual ", id,":\n", sep="")
		print(round(res,4))
		return(invisible(res))
	}
	else res
}


.likelihood_twopoint_general <- function(x, theta, dat, logbase=NULL) {

	.peel <- function(dat, sub) { 
		piv = sub[[1]]; far = sub[[2]];  mor = sub[[3]]; 
		piv_ind = .Internal(match(piv, sub[-1], nomatch=0, NULL))
		nonpiv.offs = sub[-c(1:3, piv_ind+1)] #offspring except possible pivot 
		
		farh = dat[[c(far, 1)]]; farp=dat[[c(far, 2)]]; morh=dat[[c(mor, 1)]]; morp=dat[[c(mor, 2)]]
		fa_len = length(farp); mo_len = length(morp);		
		likel = farp %*% .Internal(t.default( morp ))

		for (b in nonpiv.offs) {
			mm = seq_along(likel); dim(mm) = dim(likel)
			for(i in seq_along(farp)) for (j in seq_along(morp)) {
				mm[i,j] = (.transmitGENERAL(farh[, i], dat[[c(b,1)]][c(1,3),,drop=F], theta) * .transmitGENERAL(morh[, j], dat[[c(b,1)]][c(2,4),,drop=F], theta)) %*% dat[[c(b,2)]]
			}
			likel <-  likel * mm
		}
			
		if (piv_ind==0) return(sum(likel))	# || all(likel==0))	
		else if (piv_ind==1) res = .Internal(rowSums(likel, length(farp), length(morp), FALSE))
		else if (piv_ind==2) res = .Internal(colSums(likel, length(farp), length(morp), FALSE))
		else {
			pivh = dat[[c(piv, 1)]]; pivp = dat[[c(piv, 2)]]; pi_len = length(pivp)
			T = numeric(fa_len * mo_len * pi_len);	dim(T) = c(fa_len, mo_len, pi_len)
			for(i in seq_len(fa_len)) for(j in seq_len(mo_len))
				T[i,j,] <- .transmitGENERAL(farh[, i], pivh[c(1,3), ,drop=F], theta) * .transmitGENERAL(morh[, j], pivh[c(2,4), ,drop=F], theta)
				
			arr = .Internal(as.vector(T, 'any')) * .Internal(as.vector(likel, 'any')); dim(arr) = dim(T)
			res = .Internal(colSums(arr, fa_len*mo_len, pi_len, FALSE)) #sum for each entry of haps[[piv]]
			res = res * pivp
		}
		dat[[piv]] = list(hap = dat[[c(piv, 1)]][, res>0, drop=F], prob = res[res>0])
		return(dat)
	}
	
	if(x$hasLoops) stop("Unbroken loops in pedigree.")
	#nInd=x[['nInd']]; ped=x[['pedigree']]; chr=x[['model']][['chrom']]
	
	if(is.null(dups <- x$loop_breakers)) {
		for (sub in x$subnucs) 	{
			if (any(unlist(lapply(dat, function(dati) length(dati$prob) == 0)))) return(ifelse(is.numeric(logbase), -Inf, 0))
			dat = .peel(dat, sub)
		}
		likelihood = dat
	}
	else {
		two2one = function (matr) 1000 * matr[1, ] + matr[2, ]
		origs = match(dups[, 1], x$orig.ids); copies = match(dups[, 2], x$orig.ids)

		#For each orig, find the indices of its haplos (in orig$hap) that also occur in its copy. Then take cross product of these vectors.
		sumover = .my.grid(lapply(seq_along(origs), function(i) which(two2one(dat[[c(origs[i], 1)]]) %in% two2one(dat[[c(copies[i], 1)]]))), as.list=TRUE)
		likelihood = 0
		for(r in sumover) { #r a vector of indices: r[i] gives a column number of the hap matrix of orig[i].
			dat1 = dat; 
			for(i in seq_along(origs)) {
				origdat = dat[[origs[i]]]; 
				hap = origdat[[1]][, r[i], drop=F]
				dat1[[origs[i]]] = list(hap = hap, prob = origdat[[2]][r[i]])
				dat1[[copies[i]]] = list(hap = hap, prob = 1)
			}
			for (sub in x$subnucs) 	{
				dat1 = .peel(dat1, sub)
				if (is.numeric(dat1)) {likelihood = likelihood + dat1; break}
				if (any(unlist(lapply(dat1, function(dat1i) length(dat1i$prob) == 0)))) break #if impossible data - break out of ES-algorithm and go to next r in sumover.
			}
		}
	}
	if (is.numeric(logbase)) log(likelihood, logbase) else likelihood
}


.twopoint_general_startdata = function(x, marker1, marker2) {
	startprob <- function(h, afreq1, afreq2, founder) {
		prob = rep.int(1, ncol(h))
		if (founder) {
			m1_1 = h[1,]; m1_2 = h[2,]; m2_1 = h[3,]; m2_2 = h[4,]
			prob = prob * afreq1[m1_1] * afreq1[m1_2] * afreq2[m2_1] * afreq2[m2_2] * (((m1_1 != m1_2) | (m2_1 != m2_2)) + 1) #multiply with two if heteroz for at least 1 marker. If heteroz for both, then both phases are included in h, hence the factor 2 (not 4) in this case as well.
		}
		as.numeric(prob)
  	}
	
	## MAIN ###
	marker1 = .reduce_alleles(marker1); marker2 = .reduce_alleles(marker2)
	m1_list = .eliminate(x, marker1); 	m2_list = .eliminate(x, marker2)
	
	if (any(unlist(lapply(c(m1_list, m2_list), length)) == 0)) return(rep(list(list(prob=numeric(0))), x$nInd))	
	dat = lapply(1:x$nInd, function(i) {
		gl1 = ncol(m1_list[[i]]); gl2 = ncol(m2_list[[i]])
		hap = rbind(m1_list[[i]][, rep(1:gl1, each=gl2), drop=F ], m2_list[[i]][, rep(1:gl2, times=gl1), drop=FALSE]) #matrix with four rows: m1_1, m1_2, m2_1, m2_2
		if(i %in% x$founders) {  #doubly heterozygous founders: Include the other phase as well. (This is necessary since .eliminate returns unordered genotypes for founders.)
			doublyhet = hap[, hap[1,]!=hap[2,] & hap[3,]!=hap[4,], drop=FALSE]
			hap = cbind(hap, doublyhet[c(1,2,4,3),]) 
		}
		prob = startprob(hap, afreq1=attr(marker1, 'afreq'), afreq2=attr(marker2, 'afreq'), founder=(i %in% x$founders))
		list(hap = hap[, prob > 0, drop=F], prob = as.numeric(prob[prob > 0]))
	})
	dat
}

.transmitGENERAL <- function(parent.haps, gamete.hap, theta) { # parent.haps = c(M1_1, M1_2, M2_1, M2_2)
	norechap1 = parent.haps[c(1,3)]
	norechap2 = parent.haps[c(2,4)]
	rechap1 = parent.haps[c(1,4)]
	rechap2 = parent.haps[c(2,3)]
	if(is.matrix(gamete.hap)) 
		unlist(lapply(seq_len(ncol(gamete.hap)), function(kol) {
			gamhap = gamete.hap[,kol]
			sum(c(all(norechap1 == gamhap), all(norechap2 == gamhap), all(rechap1 == gamhap), all(rechap2 == gamhap)) 
			* c(1-theta, 1-theta, theta, theta))/2
		}))
	else {
		gamhap = gamete.hap
		sum(c(all(norechap1 == gamhap), all(norechap2 == gamhap), all(rechap1 == gamhap), all(rechap2 == gamhap)) 
			* c(1-theta, 1-theta, theta, theta))/2
	}
}


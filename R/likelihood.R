likelihood = function(x, ...)  UseMethod("likelihood", x)

likelihood.singleton = function(x, locus1, logbase=NULL, ...) {
   if(is.null(locus1) || all(locus1==0)) 
      return(if (is.numeric(logbase)) 0 else 1)
   
   m = locus1
   chrom = as.integer(attr(m, 'chrom'))
   afreq = attr(m, 'afreq') 
   if(identical(chrom, 23L) && x$pedigree[, 'SEX']==1) {# X chrom and male
      if(all(m>0) && m[1]!=m[2]) stop("Heterozygous genotype detected for X-linked marker in male individual.")
      res = afreq[m[1]]
   }
   else 
      if (0 %in% m) 
         res = afreq[m[m!=0]]
      else 
         res = prod(afreq[m]) * ifelse(m[1]!=m[2], 2, 1)  # assumes HWE
   return(if (is.numeric(logbase)) log(res, logbase) else res)
}


likelihood.list = function(x, locus1, ...) {
   if (all(sapply(x, inherits, what=c('singleton', 'linkdat'))))
      prod(sapply(1:length(x), function(i) likelihood(x[[i]], locus1[[i]], ...)))
   else stop("First argument must be either a 'linkdat' object, a 'singleton' object, or a list of such")
}

likelihood.linkdat <- function(x, locus1, locus2=NULL, theta=NULL, startdata=NULL, eliminate=0, logbase=NULL, ...) {
	
	if(x$hasLoops) stop("Unbroken loops in pedigree.")
	
	inform_subnucs = x$subnucs
	chrom = if(identical(attr(locus1, 'chrom'), 23)) 'X' else 'AUTOSOMAL'
	SEX = x$pedigree[, 'SEX']
	
	if(is.null(locus2)) {
		if(is.null(startdata)) {
			inform = .informative(x, locus1)
			inform_subnucs = inform$subnucs 
			x$founders = c(x$founders, inform$newfounders); 
			x$nonfounders = .mysetdiff(x$nonfounders, inform$newfounders)
			dat = .startdata_M(x, marker=locus1, eliminate=eliminate)
		}
		else dat = startdata
		peelFUN = function(dat, sub) .peel_M(dat, sub, chrom, SEX)	
	}
	else if(is.character(locus2)) {
		dat = if(is.null(startdata)) .startdata_MD(x, marker=locus1, eliminate=eliminate) else startdata
		peelFUN = function(dat, sub) .peel_MD(dat, sub, theta, chrom, SEX)
	}
	else {
		if(is.null(startdata)) {
			inform = .informative(x, locus1, locus2); 
			inform_subnucs = inform$subnucs; 
			x$founders = c(x$founders, inform$newfounders); 
			x$nonfounders = .mysetdiff(x$nonfounders, inform$newfounders)
			dat = .startdata_MM(x, marker1=locus1, marker2=locus2, eliminate=eliminate)
		}
		else dat = startdata
		peelFUN = function(dat, sub) .peel_MM(dat, sub, theta)
	}
	
	if (attr(dat, 'impossible')) return(ifelse(is.numeric(logbase), -Inf, 0))
	
	if(is.null(dups <- x$loop_breakers)) {
		for (sub in inform_subnucs) 	{
			dat = peelFUN(dat, sub)
			if (sub$pivtype > 0 && attr(dat, 'impossible')) return(ifelse(is.numeric(logbase), -Inf, 0))
		}
		likelihood = dat
	}
	else {
		two2one = function(matr) if(is.matrix(matr)) 1000 * matr[1, ] + matr[2, ] else matr  #if input is vector (i.e. X-linked male genotypes), return it unchanged
		origs = match(dups[, 1], x$orig.ids); copies = match(dups[, 2], x$orig.ids)

		#For each orig, find the indices of its haplos (in orig$hap) that also occur in its copy. Then take cross product of these vectors.
		loopgrid = .my.grid(lapply(seq_along(origs), function(i)  { ori = two2one(dat[[c(origs[i], 1)]]); seq_along(ori)[ori %in% two2one(dat[[c(copies[i], 1)]])] }), as.list=TRUE)
		likelihood = 0
		for(r in loopgrid) { #r a vector of indices: r[i] gives a column number of the hap matrix of orig[i].
			dat1 = dat; attr(dat1, 'impossible') = FALSE
			for(i in seq_along(origs)) {
				orig.int = origs[i]; copy.int = copies[i]
				origdat = dat[[orig.int]]; 
				hap = if(chrom=='X' && SEX[orig.int]==1) origdat[[1]][r[i]] else origdat[[1]][, r[i], drop=F]
				prob = origdat[[2]][r[i]]; if(sum(prob)==0) print("Loop-loekke: Alle sannsynligheter er null. Magnus lurer paa om dette kan gi feilmelding.")
				dat1[[orig.int]] = list(hap = hap, prob = prob)
				dat1[[copy.int]] = list(hap = hap, prob = 1)
			}
			
			for (sub in inform_subnucs) 	{
				pivtype = sub$pivtype
				dat1 = peelFUN(dat1, sub)
				if (pivtype > 0 && attr(dat1, 'impossible')) {break} #if impossible data - break out of ES-algorithm and go to next r in loopgrid.
				if (pivtype == 0) likelihood = likelihood + dat1
			}
		}
	}
	if (is.numeric(logbase)) log(likelihood, logbase) else likelihood
}

#### FUNCTIONS FOR CREATING THE INTITIAL HAPLOTYPE COMBINATIONS W/PROBABILITIES.

.startdata_M = function(x, marker, eliminate=0) {
	marker = .reduce_alleles(marker)

	afreq = attr(marker, 'afreq')
	chromX = identical(attr(marker, 'chrom'), 23)

	impossible = FALSE
	if(chromX) {
		glist = .build_genolist_X(x, marker, eliminate)
		if (attr(glist, 'impossible')) {dat=list(); attr(dat, 'impossible')=TRUE; return(dat)}
		sex = x$pedigree[, 'SEX']
		dat = lapply(1:x$nInd, function(i) {
			h = glist[[i]]
			if(i %in% x$founders) {
				prob = switch(sex[i], afreq[h], afreq[h[1,]] * afreq[h[2,]] * ((h[1,] != h[2,]) + 1))
				if(sum(prob)==0) impossible = TRUE
			}
			else  prob = rep.int(1, length(h)/sex[i])	
			list(hap = h, prob = as.numeric(prob))
		})
	}
	else {
		glist = .build_genolist(x, marker, eliminate)
		if (attr(glist, 'impossible')) {dat=list(); attr(dat, 'impossible')=TRUE; return(dat)}
		dat = lapply(1:x$nInd, function(i) {
			h = glist[[i]]
			if(i %in% x$founders) {
				prob = afreq[h[1,]] * afreq[h[2,]] * ((h[1,] != h[2,]) + 1)
				if(sum(prob)==0) impossible = TRUE
			}
			else  prob = rep.int(1, ncol(h))	
			list(hap = h, prob = as.numeric(prob))
		})
	}
	attr(dat, 'impossible') = impossible
	dat
}


.startdata_MD = function(x, marker, eliminate=0) {
	startprob <- function(h, model, afreq, aff, founder) {
 		pat = h[1,]; mat = h[2,]
		al1 = abs(pat); al2 = abs(mat)
		d.no = (pat<0) + (mat<0)
		prob = switch(aff+1, rep.int(1, length(d.no)), (1 - model$penetrances)[d.no+1], model$penetrances[d.no+1])
		if (founder)
			prob = prob * afreq[al1] * afreq[al2] * ((al1 != al2) + 1) * model$dfreq^d.no * (1-model$dfreq)^(2-d.no)
		as.numeric(prob)
  	}
	
	startprob_X <- function(h, model, afreq, sex, aff, founder) {
 		switch(sex, {mat = h #vector
			d.no = as.numeric(mat<0)
			prob <- switch(aff+1, rep.int(1, length(d.no)), (1 - model$penetrances$male)[d.no+1], model$penetrances$male[d.no+1])
			if (founder)  
				prob <- prob * afreq[abs(mat)] * c(1-model$dfreq, model$dfreq)[d.no+1]
		}, {pat = h[1,]; mat = h[2,]
			al1 = abs(pat); al2 = abs(mat)
			d.no = (pat<0) + (mat<0)
			prob <- switch(aff+1, rep.int(1, length(d.no)), (1 - model$penetrances$female)[d.no+1], model$penetrances$female[d.no+1])
			if (founder)
				prob <- prob * afreq[al1] * afreq[al2] * ((al1 != al2) + 1) * model$dfreq^d.no * (1-model$dfreq)^(2-d.no)
		})
		as.numeric(prob)
  	}
	
	## MAIN ###
	marker = .reduce_alleles(marker)
	afreq = attr(marker, 'afreq')
	chromX = identical(attr(marker, 'chrom'), 23)
	AFF = x$pedigree[, 'AFF']; FOU = (1:x$nInd) %in% x$founders
	dlist = cbind(c(-1,-1),c(-1,1),c(1,-1),c(1,1)) #D=-1; N=1
	impossible = FALSE
	if(chromX) {
		glist = .build_genolist_X(x, marker, eliminate)
		if (attr(glist, 'impossible')) {dat=list(); attr(dat, 'impossible')=TRUE; return(dat)}

		SEX = x$pedigree[,'SEX']
		dat = lapply(1:x$nInd, function(i) {
			switch(SEX[i], 
			{ hap = c(glist[[i]], -glist[[i]])
			  prob = startprob_X(hap, model=x$model, afreq=attr(marker, 'afreq'), sex=1, aff=AFF[i], founder=FOU[i])
			  list(hap = hap[prob > 0], prob = prob[prob > 0])},
			{ gl = ncol(glist[[i]])  
			  hap = glist[[i]][, rep(1:gl, each=4), drop=F ] *  dlist[, rep(1:4, times=gl), drop=F]
			  prob = startprob_X(hap, model=x$model, afreq=attr(marker, 'afreq'), sex=2, aff=AFF[i], founder=FOU[i])
			  if (all(prob==0)) impossible = TRUE
			  list(hap = hap[, prob > 0, drop=F], prob = as.numeric(prob[prob > 0]))
			})
		})	
	}
	else {
		glist = .build_genolist(x, marker, eliminate)
		if (attr(glist, 'impossible')) {dat=list(); attr(dat, 'impossible')=TRUE; return(dat)}
	
		dat = lapply(1:x$nInd, function(i) {
			gl = ncol(glist[[i]])  
			hap = glist[[i]][, rep(1:gl, each=4), drop=F ] *  dlist[, rep(1:4, times=gl), drop=F ]
			prob = startprob(hap, model=x$model, afreq=attr(marker, 'afreq'), aff=AFF[i], founder=FOU[i])
			if(sum(prob)==0) impossible=TRUE
			list(hap = hap[, prob > 0, drop=F], prob = as.numeric(prob[prob > 0]))
		})
	}
	attr(dat, 'impossible') = impossible
	dat
}



.startdata_MM = function(x, marker1, marker2, eliminate=0) {
	startprob <- function(h, afreq1, afreq2, founder) {
		prob = rep.int(1, ncol(h))
		if (founder) {
			m1_1 = h[1,]; m1_2 = h[2,]; m2_1 = h[3,]; m2_2 = h[4,]
			prob = prob * afreq1[m1_1] * afreq1[m1_2] * afreq2[m2_1] * afreq2[m2_2] * (((m1_1 != m1_2) | (m2_1 != m2_2)) + 1) #multiply with two if heteroz for at least 1 marker. If heteroz for both, then both phases are included in h, hence the factor 2 (not 4) in this case as well.
		}
		as.numeric(prob)
  	}

	marker1 = .reduce_alleles(marker1); marker2 = .reduce_alleles(marker2)
	m1_list = .build_genolist(x, marker1, eliminate)
 	m2_list = .build_genolist(x, marker2, eliminate)
	impossible = FALSE
	
	if (any(unlist(lapply(c(m1_list, m2_list), length)) == 0)) return(rep(list(list(prob=numeric(0))), x$nInd))	
	dat = lapply(1:x$nInd, function(i) {
		gl1 = ncol(m1_list[[i]]); gl2 = ncol(m2_list[[i]])
		hap = rbind(m1_list[[i]][, rep(1:gl1, each=gl2), drop=F ], m2_list[[i]][, rep(1:gl2, times=gl1), drop=FALSE]) #matrix with four rows: m1_1, m1_2, m2_1, m2_2
		if(i %in% x$founders) {  #doubly heterozygous founders: Include the other phase as well. (This is necessary since .build_genolist returns unordered genotypes for founders.)
			doublyhet = hap[, hap[1,]!=hap[2,] & hap[3,]!=hap[4,], drop=FALSE]
			hap = cbind(hap, doublyhet[c(1,2,4,3),]) 
		}
		prob = startprob(hap, afreq1=attr(marker1, 'afreq'), afreq2=attr(marker2, 'afreq'), founder=(i %in% x$founders))
		if(sum(prob)==0) impossible=TRUE
		list(hap = hap[, prob > 0, drop=F], prob = as.numeric(prob[prob > 0]))
	})
	attr(dat, 'impossible') = impossible
	dat
}




#### PEELING FUNCTIONS (one for each case: single Marker, Marker-Disease, Marker-Marker)


.peel_M <- function(dat, sub, chrom, SEX) { 
	far = sub[['father']]; mor = sub[['mother']]; offs = sub[['offspring']]; piv = sub[['pivot']]; pivtype = sub[['pivtype']]
	if(pivtype==3) offs = offs[offs != piv] #pivtype indicates who is pivot: 0 = none; 1 = father; 2 = mother; 3 = an offspring
	farh = dat[[c(far, 1)]]; morh = dat[[c(mor, 1)]]
	likel = dat[[c(far, 2)]] %*% t.default(dat[[c(mor, 2)]])
	dims = dim(likel);	fa_len = dims[1L];  mo_len = dims[2L]		
	
	.trans_M <- function(parent, childhap) unlist(lapply(seq_len(ncol(parent)), function(i) ((parent[1,i] == childhap) + (parent[2,i] == childhap))/2)) #parent = matrix with 2 rows; childhap = vector of any length (parental allele)

	switch(chrom, 
	AUTOSOMAL = {
		for (datb in dat[offs]) {
			bh = datb[[1]]; bp = datb[[2]]; bl = length(bp)
			trans_pats = .trans_M(farh, bh[1, ])
			trans_mats = .trans_M(morh, bh[2, ])
			dim(trans_mats) = c(bl, mo_len); trans_mats_rep = as.numeric(do.call(rbind, rep(list(trans_mats), fa_len)))
			mm = .colSums((trans_pats * bp) * trans_mats_rep, bl, fa_len*mo_len)
			likel = likel * mm
		}
		if(pivtype==0) return(sum(likel))
		
		res = switch(pivtype, .rowSums(likel, fa_len, mo_len), .colSums(likel, fa_len, mo_len),
			{	pivh = dat[[c(piv, 1)]]; pivp = dat[[c(piv, 2)]]; pi_len = length(pivp)
				T = numeric(fa_len * mo_len * pi_len);	dim(T) = c(fa_len, mo_len, pi_len)
				trans_pats = .trans_M(farh, pivh[1, ]); dim(trans_pats) = c(pi_len, fa_len);
				trans_mats = .trans_M(morh, pivh[2, ]); dim(trans_mats) = c(pi_len, mo_len);
				for(i in seq_len(fa_len)) {
					transpat = trans_pats[, i]
					for(j in seq_len(mo_len)) T[i,j,] = transpat * trans_mats[, j]
				}
				arr = as.vector(T) * as.vector(likel); dim(arr) = dim(T)
				res = .colSums(arr, fa_len*mo_len, pi_len) #sum for each entry of haps[[piv]]
				res * pivp
		  })
		pivhap_update = dat[[c(piv, 1)]][, res>0, drop=F]
	},
	X = {
		for (b in offs) {
			datb = dat[[b]]; bh = datb[[1]]; bp = datb[[2]]; bl = length(bp)
			switch(SEX[b],
			{ trans_mats = .trans_M(morh, bh)
			  mm = rep(.colSums(trans_mats * bp, bl, mo_len), each=fa_len); 
			},
			{ trans_pats = unlist(lapply(farh, function(fh) as.numeric(fh == bh[1, ])))
			  trans_mats = .trans_M(morh, bh[2, ])
			  dim(trans_mats) = c(bl, mo_len); trans_mats_rep = as.numeric(do.call(rbind, rep(list(trans_mats), fa_len))) #TODO: improve with 'rep'?
			  mm = .colSums((trans_pats * bp) * trans_mats_rep, bl, fa_len*mo_len)
			})
			likel = likel * mm 
		}
		if(pivtype==0) return(sum(likel))
			
		res = switch(pivtype, .rowSums(likel, fa_len, mo_len), .colSums(likel, fa_len, mo_len),
			{	pivh = dat[[c(piv, 1)]]; pivp = dat[[c(piv, 2)]]; pi_len = length(pivp)
				switch(SEX[piv],
				{	trans_mats = .trans_M(morh, pivh)
					T = rep(trans_mats, each=fa_len)
				}, 
				{	T = numeric(fa_len * mo_len * pi_len);	dim(T) = c(fa_len, mo_len, pi_len)
					trans_mats = .trans_M(morh, pivh[2, ]); dim(trans_mats) = c(pi_len, mo_len)
					for (i in seq_len(fa_len))
						T[i,,] = t.default(as.numeric(farh[i] == pivh[1,]) * trans_mats) #TODO:make faster?
				})
				arr = as.vector(T) * as.vector(likel); dim(arr) = dim(T)
				res = .colSums(arr, fa_len*mo_len, pi_len) #sum for each entry of haps[[piv]]
				res * pivp
			})
		pivhap_update = switch(SEX[piv], dat[[c(piv, 1)]][res>0], dat[[c(piv, 1)]][, res>0, drop=F])
	})
	dat[[piv]] = list(hap = pivhap_update, prob = res[res>0])
	if(sum(res)==0) attr(dat, 'impossible') = TRUE
	return(dat)
}

.peel_MD <- function(dat, sub, theta, chrom, SEX) {
	far = sub[['father']]; mor = sub[['mother']]; offs = nonpiv.offs = sub[['offspring']]; piv = sub[['pivot']]; pivtype = sub[['pivtype']]
	if(pivtype==3) non.pivoffs = offs[offs != piv]
	farh = dat[[c(far, 1)]]; morh = dat[[c(mor, 1)]]
	likel = dat[[c(far, 2)]] %*% t.default(dat[[c(mor, 2)]])
	dims = dim(likel);	fa_len = dims[1L];  mo_len = dims[2L]		
	
	.trans_MD <- function(parent_ph, child_haplo, theta) {
		if (theta==0) 	return(((parent_ph[1] == child_haplo) + (parent_ph[2] == child_haplo))/2)
		parent_rec = abs(parent_ph) * sign(parent_ph[2:1])
		((parent_ph[1]==child_haplo) * (1-theta) + (parent_ph[2]==child_haplo) * (1-theta) 	+  (parent_rec[1]==child_haplo) * theta + (parent_rec[2]==child_haplo) * theta)/2
	}

	switch(chrom, AUTOSOMAL = {
		alloffs_pat = unique.default(unlist(lapply(offs, function(i) dat[[c(i, 1)]][1,])))
		alloffs_mat = unique.default(unlist(lapply(offs, function(i) dat[[c(i, 1)]][2,])))
		fathertrans = unlist(lapply(seq_len(fa_len), function(i) .trans_MD(farh[, i], alloffs_pat, theta)))
		mothertrans = unlist(lapply(seq_len(mo_len), function(j) .trans_MD(morh[, j], alloffs_mat, theta)))
		dim(fathertrans) = c(length(alloffs_pat), fa_len); dim(mothertrans) = c(length(alloffs_mat), mo_len)
		mm_init = numeric(fa_len*mo_len); dim(mm_init) = dims
		for (b in nonpiv.offs) {
			mm = mm_init; 
			bh = dat[[c(b,1)]]; bp = dat[[c(b,2)]]
			b_pat = match(bh[1, ], alloffs_pat);	b_mat = match(bh[2, ], alloffs_mat)
			for(i in seq_len(fa_len))  for (j in seq_len(mo_len))
					mm[i,j] = (fathertrans[b_pat, i] * mothertrans[b_mat, j]) %*% bp #TODO: make faster
			likel = likel * mm
		}
		if(pivtype==0) return(sum(likel))
		
		res = switch(pivtype, .rowSums(likel, fa_len, mo_len), .colSums(likel, fa_len, mo_len),
		{
			pivh = dat[[c(piv, 1)]]; pivp = dat[[c(piv, 2)]]; pi_len = length(pivp)
			piv_pat = match(pivh[1, ], alloffs_pat);	piv_mat = match(pivh[2, ], alloffs_mat)
			T = numeric(fa_len * mo_len * pi_len);	dim(T) = c(fa_len, mo_len, pi_len)
			for(i in seq_len(fa_len)) for(j in seq_len(mo_len)) #TODO: splitt up to make faster!
				T[i,j,] <- fathertrans[piv_pat, i] * mothertrans[piv_mat, j]	
			arr = as.vector(T) * as.vector(likel); dim(arr) = dim(T)
			res = .colSums(arr, fa_len*mo_len, pi_len) #sum for each entry of haps[[piv]]
			res = res * pivp
		})
		pivhap_update = dat[[c(piv, 1)]][, res>0, drop=F]
	},
	X = {
		mm_init = numeric(length(likel)); dim(mm_init) = dims
		for (b in nonpiv.offs) {
			mm = mm_init
			switch(SEX[b], for (j in seq_len(mo_len))  mm[, j] = .trans_MD(morh[, j], dat[[c(b,1)]], theta) %*% dat[[c(b,2)]], 
			{for (j in seq_len(mo_len)) {
					transmother = .trans_MD(morh[, j], dat[[c(b,1)]][2,], theta)
					for (i in seq_len(fa_len))
						mm[i,j] = (as.numeric(farh[i] == dat[[c(b,1)]][1,]) * transmother) %*% dat[[c(b,2)]]	
			}})
			likel = likel * mm
		}
		if(pivtype==0) return(sum(likel))
		
		res = switch(pivtype, .rowSums(likel, fa_len, mo_len), .colSums(likel, fa_len, mo_len),
		{
			pivh = dat[[c(piv, 1)]]; pivp = dat[[c(piv, 2)]]; pi_len = length(pivp)
			T = numeric(fa_len * mo_len * pi_len);	dim(T) = c(fa_len, mo_len, pi_len)
			switch(SEX[piv],
			{for (j in seq_len(mo_len))
				transmother = .trans_MD(morh[, j], pivh, theta)
				for(i in seq_len(fa_len)) T[i,j,] = transmother
			}, 
			{for (j in seq_len(mo_len)) {
				transmother = .trans_MD(morh[, j], pivh[2,], theta)
				for (i in seq_len(fa_len))
					T[i,j,] = (as.numeric(farh[i] == pivh[1,]) * transmother)
			}})

			arr = as.vector(T) * as.vector(likel); dim(arr) = dim(T)
			res = .colSums(arr, fa_len*mo_len, pi_len) #sum for each entry of haps[[piv]]
			res * pivp
		})
		pivhap_update = switch(SEX[piv], dat[[c(piv, 1)]][res>0], dat[[c(piv, 1)]][, res>0, drop=F])
	})
	dat[[piv]] = list(hap = pivhap_update, prob = res[res>0])
	if(all(res==0)) attr(dat, 'impossible') = TRUE
	return(dat)
}



.peel_MM <- function(dat, sub, theta) { 
	far = sub[['father']]; mor = sub[['mother']]; offs = sub[['offspring']]; piv = sub[['pivot']]; pivtype = sub[['pivtype']]
	if(pivtype==3) offs = offs[offs != piv]
	
	farh = dat[[c(far, 1)]]; morh=dat[[c(mor, 1)]]
	likel = dat[[c(far, 2)]] %*% t.default(dat[[c(mor, 2)]])
	dims = dim(likel);	fa_len = dims[1L];  mo_len = dims[2L];	
	mm_init = numeric(fa_len*mo_len); dim(mm_init) = dims

	.trans_MM <- function(parent.haps, gamete.hap, theta) { # parent.haps = c(M1_1, M1_2, M2_1, M2_2)
		norechap1 = parent.haps[c(1,3)];	norechap2 = parent.haps[c(2,4)]
		rechap1 = parent.haps[c(1,4)];  	rechap2 = parent.haps[c(2,3)]
		if(is.matrix(gamete.hap)) 
			unlist(lapply(seq_len(ncol(gamete.hap)), function(kol) {
				gamhap = gamete.hap[,kol]
				sum(c(all(norechap1 == gamhap), all(norechap2 == gamhap), all(rechap1 == gamhap), all(rechap2 == gamhap)) * c(1-theta, 1-theta, theta, theta))/2
			}))
		else {
			gamhap = gamete.hap
			sum(c(all(norechap1 == gamhap), all(norechap2 == gamhap), all(rechap1 == gamhap), all(rechap2 == gamhap)) * c(1-theta, 1-theta, theta, theta))/2
		}
	}

	for (b in offs) {
		mm = mm_init
		bh = dat[[c(b,1)]]; bp = dat[[c(b,2)]]
		for(i in seq_len(fa_len)) {
			transfather = .trans_MM(farh[, i], bh[c(1,3),,drop=F], theta)
			for (j in seq_len(mo_len))
				mm[i,j] = (transfather * .trans_MM(morh[, j], bh[c(2,4),,drop=F], theta)) %*% bp
		}
		likel <-  likel * mm
	}
	if(pivtype==0) return(sum(likel))
	
	res = switch(pivtype, .rowSums(likel, fa_len, mo_len), .colSums(likel, fa_len, mo_len),
	{ pivh = dat[[c(piv, 1)]]; pivp = dat[[c(piv, 2)]]; pi_len = length(pivp)
		T = numeric(fa_len * mo_len * pi_len);	dim(T) = c(fa_len, mo_len, pi_len)
		for(i in seq_len(fa_len)) {
			transfather = .trans_MM(farh[, i], pivh[c(1,3), ,drop=F], theta)
			for(j in seq_len(mo_len))
				T[i,j,] <- transfather * .trans_MM(morh[, j], pivh[c(2,4), ,drop=F], theta)
		}
		arr = as.vector(T) * as.vector(likel); dim(arr) = dim(T)
		res = .colSums(arr, fa_len*mo_len, pi_len) #sum for each entry of haps[[piv]]
		res = res * pivp
	})
	dat[[piv]] = list(hap = dat[[c(piv, 1)]][, res>0, drop=F], prob = res[res>0])
	return(dat)
}

###### OTHER AUXILIARY FUNCTIONS



.dropIndivs = function(x, marker, marker2=NULL, peel=FALSE, internal.ids=TRUE) { #NB: peelingen virker ikke for komplekse peds...   #outputs original IDs of indivs that can be skipped in the peeling
	if(!is.null(marker2)) marker = marker + marker2
	if(all(marker[,1] > 0)) return(numeric(0))
	miss = seq_len(x$nInd)[marker[,1]==0 & marker[,2]==0]
	p = x$pedigree
	leaves = .mysetdiff(p[, "ID"], p[, c("FID", "MID")])
	throw = .myintersect(leaves, miss)
	if(peel) {stop("peel = T is not implemented yet") #TODO fix this for complex peds...
		subs = x$subnucs
		while (length(subs)>1) {
			alloffs.uninf = unlist(lapply(subs, function(s) all(s[-c(1:3)] %in% throw)))
			if (!any(alloffs.uninf)) break
			throw = c(throw, unlist(lapply(subs[alloffs.uninf], function(s) switch(sum(parmis <- (s[2:3] %in% miss)) + 1, NULL, {mis = s[2:3][parmis]; if(s[1]!=mis) mis}, s[2:3]))))
			subs = subs[!alloffs.uninf]
		}
	}
	if(!is.null(lb <- x$loop_breakers)) throw = .mysetdiff(throw, .internalID(x, lb))
	if(internal.ids) as.numeric(throw) else as.numeric(x$orig.ids[throw])
}

.informative = function(x, marker, marker2=NULL) {
	if(!is.null(marker2)) marker = marker + marker2
	if(all(marker[,1] > 0)) return(list(subnucs=x$subnucs, newfounders=numeric(0)))
	newfounders = numeric(0)
	new_subnucs = list()
	p = x$pedigree
	is_miss = marker[,1]==0 & marker[,2]==0
	is_miss[loop_int <- .internalID(x, x$loop_breakers)] = F # works (and quick) also if no loops.
	is_uninf_leaf = is_miss & !p[, "ID"] %in% p[, c("FID", "MID")]
	is_uninf_fou = is_miss & seq_len(x$nInd) %in% x$founders
	
	for(sub in x$subnucs) {
		fa = sub[['father']]; mo = sub[['mother']]; offs = sub[['offspring']]; pivot = sub[['pivot']]
		sub[['offspring']] = offs[!is_uninf_leaf[offs]]
		noffs = length(sub[['offspring']])
		switch(sub[['pivtype']], 
			{if(noffs==0 && is_uninf_fou[mo]) sub = NULL}, 
			{if(noffs==0 && is_uninf_fou[fa]) sub = NULL},
			{if(noffs==1 && is_uninf_fou[fa] && is_uninf_fou[mo]) {newfounders = c(newfounders, pivot); sub = NULL}})
		if(!is.null(sub)) {
         new_subnucs = c(new_subnucs, list(sub))
         is_uninf_fou[pivot] = FALSE  #added in v0.8-1 to correct a bug marking certain "middle" subnucs uninformative
      }
   }
	list(subnucs = new_subnucs, newfounders = newfounders)
}


.make.grid.subset = function(x, partialmarker, ids, chrom, make.grid=T) {
	int.ids = .internalID(x, ids)
	nall = attr(partialmarker, 'nalleles')
   if(missing(chrom)) 
      chrom = if(identical(attr(partialmarker, 'chrom'), 23)) 'X' else 'AUTOSOMAL'

	allg = .allGenotypes(nall)
	allg_ref = 1000*(allg[,1] + allg[,2]) + abs(allg[,1] - allg[,2])
	
   match_ref_rows = function(genomatr) # In: matrix with 2 rows (each column a genotype). Out: vector of 'allg' row numbers
      sort.int(unique.default(match(1000*(genomatr[1,] + genomatr[2,]) + abs(genomatr[1,] - genomatr[2,]), allg_ref)))
   
   switch(chrom, 
	'AUTOSOMAL' = {
		glist = .build_genolist(x, partialmarker, eliminate=100)
		if (attr(glist, 'impossible')) stop("Impossible partial marker")
		rows = lapply(glist[int.ids], match_ref_rows)
	},
	'X' = {
		SEX = x$pedigree[, 'SEX']
		glist = .build_genolist_X(x, partialmarker, eliminate=100)
		if (attr(glist, 'impossible')) stop("Impossible partial marker")
		rows = lapply(int.ids, function(i) switch(SEX[i], glist[[i]], match_ref_rows(glist[[i]])))
	})
	if (make.grid) .my.grid(rows) else rows
}

				
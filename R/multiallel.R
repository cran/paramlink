
.trans <- function(parent_ph, child_haplo, theta) {
	if (theta==0)
		return(((parent_ph[1] == child_haplo) + (parent_ph[2] == child_haplo))/2)
	parent_rec = abs(parent_ph) * sign(parent_ph[2:1])

	#sum((c(parent_ph, parent_rec) == child_haplo) * c(1-theta, 1-theta, theta, theta))/2
	((parent_ph[1]==child_haplo) * (1-theta) + (parent_ph[2]==child_haplo) * (1-theta)
	+ (parent_rec[1]==child_haplo) * theta + (parent_rec[2]==child_haplo) * theta)/2
}

.likelihood_multi <- function(x, theta, dat, logbase=NULL) {

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
				mm[i,j] = (.trans(farh[, i], dat[[c(b,1)]][1,], theta) * .trans(morh[, j], dat[[c(b,1)]][2,], theta)) %*% dat[[c(b,2)]]
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
				T[i,j,] <- .trans(farh[, i], pivh[1, ], theta) * .trans(morh[, j], pivh[2, ], theta)
				
			arr = .Internal(as.vector(T, 'any')) * .Internal(as.vector(likel, 'any')); dim(arr) = dim(T)
			res = .Internal(colSums(arr, fa_len*mo_len, pi_len, FALSE)) #sum for each entry of haps[[piv]]
			res = res * pivp
		}
		dat[[piv]] = list(hap = dat[[c(piv, 1)]][, res>0, drop=F], prob = res[res>0])
		return(dat)
	}
	
	if(x$hasLoops) stop("Unbroken loops in pedigree.")
	nInd=x[['nInd']]; ped=x[['pedigree']]; chr=x[['model']][['chrom']]
	if (any(unlist(lapply(dat, function(dati) length(dati$prob) == 0)))) return(ifelse(is.numeric(logbase), -Inf, 0))
	
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
		sumover = .my.grid(lapply(seq_along(origs), function(i)
			which(two2one(dat[[c(origs[i], 1)]]) %in% two2one(dat[[c(copies[i], 1)]]))), as.list=TRUE)
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
#one2two = function(u) rbind(v1 <- round(u/1000), u - v1*1000)

.lodm <- function(x, markers=seq_len(x$nMark), theta=0, loop_breakers=NULL, max.only=FALSE, verbose=FALSE) {

	stopifnot(class(x)=="linkdat")
	if (is.null(x$model)) 	stop("No model set.")
	if (x$nMark==0) 		stop("No marker data indicated.")
	if (x$hasLoops)		{	
		if(is.null(loop_breakers)) stop("The pedigree has loops. Please indicate loop breakers.")
		x = breakLoops(x, loop_breakers)
	}
	
	theta_list = as.list(theta)
	lods = vapply(x$markerdata[markers], function(marker) {
		start.dat = .likelihood_startdata(x, marker)
		denom = .likelihood_multi(x=x, theta=0.5, dat=start.dat, logbase=10)
		unlist(lapply(theta_list, function(theta) .likelihood_multi(x=x, theta=theta, dat=start.dat, logbase=10) - denom))
	}, FUN.VALUE=numeric(length(theta)))
	
	map = .getMap(x, na.action=1)[markers, ,drop=FALSE]

	res = structure(lods, dim=c(length(theta), length(markers)), dimnames = list(theta, map$MARKER), analysis="mlink", map=map, class="linkres")
	if (verbose) summary(res)
	if (max.only) ifelse(all(is.na(res)), NA, max(res, na.rm=TRUE))
	else res
}


.eliminate <- function(x, mm) {  #mm: marker matrix, dim = (nInd , 2)
	nallel = attr(mm, 'nalleles')
	genolist = lapply(1:x$nInd, function(i)  {
		a = mm[i,1]; b = mm[i,2]
		if(a == 0) a = seq_len(nallel)
		if(b == 0) b = seq_len(nallel)
		m = matrix( c(rep(a, each=length(b)), rep.int(b, times=2*length(a)), rep(a, each=length(b))), byrow=T, nrow=2)
		if ((i %in% x$founders) && (!x$orig.ids[i] %in% x$loop_breakers)) 
			m = m[, m[1,] <= m[2,], drop=FALSE]
		unique.matrix(m, MARGIN=2)
	})
 
	for (sub in x$subnucs) {
		far = genolist[[ sub[[2]] ]]; mor = genolist[[ sub[[3]] ]]; barnlist=sub[-c(1:3)]
		far1 = far[1,]; far2 = far[2,]; fl = length(far1)
		mor_ext1 = rep(mor[1,], each=fl); mor_ext2 = rep(mor[2,], each=fl); ml = ncol(mor);
		
		comp <- lapply( genolist[ barnlist ], function(b) {
			b1 = rep(b[1,], each=fl*ml); b2 = rep(b[2,], each=fl*ml)
			compv = (far1 == b1 | far2 == b1) & (mor_ext1 == b2 | mor_ext2 == b2) #far and mor_ext recycle
			dim(compv) = c(fl,ml,ncol(b))
			compv
		})			
		parent_keep = rep.int(T, fl*ml)
		for(arr in comp) 
			parent_keep = parent_keep * rowSums(arr, dims=2) #sums over dimension 3, i.e. parent_keep goes to 0 if there is an offsp with no compat genotypes. Faster than apply(,3,any)
		far_keep = rowSums(parent_keep)>0
		mor_keep = colSums(parent_keep)>0
		genolist[[ sub[[2]] ]] <- far[, far_keep, drop=F]
		genolist[[ sub[[3]] ]] <- mor[, mor_keep, drop=F]
		offs_keep <- lapply(comp, function(arr) colSums(arr[far_keep, mor_keep,,drop=F], dims=2)>0) #sums over dimensions 1:2
		for (bno in seq_along(barnlist)) genolist[[barnlist[bno]]] <- genolist[[barnlist[bno]]][ , offs_keep[[bno]], drop=F]
	}
	genolist
}

.reduce_alleles = function(marker) {
	if(all(marker!=0)) return(marker)
	present = .Internal(unique(as.numeric(marker),F,F))
	if(length(present) >= attr(marker, 'nalleles')) return(marker)  # returns if at most one allele is not present
	present_alleles = present[present > 0]
	present_freq = attr(marker, 'afreq')[present_alleles]
	new_marker = match(marker, present_alleles, nomatch=0)
	attributes(new_marker) = attributes(marker)
	attr(new_marker, 'alleles') = c(attr(marker, 'alleles')[present_alleles], "dummy")
	attr(new_marker, 'nalleles') = length(present)
	attr(new_marker, 'afreq') = c(present_freq, 1 - sum(present_freq))
	new_marker
}

.likelihood_startdata = function(x, marker) {
	startprob <- function(h, model, afreq, aff, founder) {
 		pat = h[1,]; mat = h[2,]
		al1 = abs(pat); al2 = abs(mat)
		d.no = (pat<0) + (mat<0)
		if (aff == 2) prob = model$penetrances[d.no+1] else prob = (1 - model$penetrances)[d.no+1]
		if (founder)
			prob = prob * afreq[al1] * afreq[al2] * ((al1 != al2) + 1) * model$dfreq^d.no * (1-model$dfreq)^(2-d.no)
		as.numeric(prob)
  	}
	
	## MAIN ###
	marker = .reduce_alleles(marker)
	glist = .eliminate(x, marker)
	if (any(unlist(lapply(glist, length)) == 0)) return(rep(list(list(prob=numeric(0))), x$nInd))
	dlist = cbind(c(-1,-1),c(-1,1),c(1,-1),c(1,1)) #D=-1; N=1
	dat = lapply(1:x$nInd, function(i) {
		gl = ncol(glist[[i]])  
		hap = glist[[i]][, rep(1:gl, each=4), drop=F ] *  dlist[, rep(1:4, times=gl), drop=F ]
		prob = startprob(hap, model=x$model, afreq=attr(marker, 'afreq'), aff=x$pedigree[i, 'AFF'], founder=(i %in% x$founder))
		list(hap = hap[, prob > 0, drop=F], prob = as.numeric(prob[prob > 0]))
	})
	dat
}



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
		attr(marker, 'chrom') = x$model$chrom
		start.dat = .startdata_MD(x, marker)
		denom = likelihood(x, locus1=marker, locus2="disease", theta=0.5, startdata = start.dat, logbase=10)
		unlist(lapply(theta_list, function(tht) likelihood(x, locus1=marker, locus2="disease", theta=tht, startdata = start.dat, logbase=10) - denom))
	}, FUN.VALUE=numeric(length(theta)))
	
	map = .getMap(x, na.action=1, verbose=F)[markers, ,drop=FALSE]

	res = structure(lods, dim=c(length(theta), length(markers)), dimnames = list(theta, map$MARKER), analysis="mlink", map=map, class="linkres")
	if (verbose) summary(res)
	if (max.only) ifelse(all(is.na(res)), NA, max(res, na.rm=TRUE))
	else res
}


.build_genolist <- function(x, marker, eliminate=0) {  #mm: marker matrix, dim = (nInd , 2)
	n = attr(marker, 'nalleles')
	nseq = seq_len(n)
	.COMPLETE = { tmp1 = rep(nseq, each=n); tmp2 = rep.int(nseq, times=n);
		fath = c(tmp1, tmp2); moth = c(tmp2, tmp1)
		rbind(fath, moth, deparse.level=0)[, !duplicated.default(fath*1000 + moth), drop=F] #faster than unique(m, MARGIN=2)
	}
	genolist = lapply(1:x$nInd, function(i)  {
		g_i = marker[i, ];
		m = switch(sum(g_i == 0) + 1, 
		{ a=g_i[1]; b=g_i[2]; if(a==b) cbind(g_i, deparse.level=0) else cbind(g_i, c(b, a), deparse.level=0)},
		{ nz = g_i[g_i!=0]; rbind(c(nseq, rep.int(nz, n-1)), c(rep.int(nz, n), nseq[-nz]), deparse.level=0)},
		{ .COMPLETE })
		if ((i %in% x$founders) && (!x$orig.ids[i] %in% x$loop_breakers)) 
			m = m[, m[1,] <= m[2,], drop=FALSE]
		m
	})
	attr(genolist, 'impossible') = FALSE
	.eliminate(x, genolist, n, repeats=eliminate)
}


.eliminate = function(x, genolist, nall, repeats=0) {
	if(repeats==0 || attr(genolist, 'impossible')) return(genolist)
	offs = lapply(1:x$nInd, function(i) offspring(x, i, original.id=FALSE))
	ncols_ny = unlist(lapply(genolist, ncol))
	p = x$pedigree
	informative = logical(x$nInd)
	for(k in seq_len(repeats)) {
		ncols = ncols_ny
		informative[x$founders] = (ncols[x$founders] < nall*(nall+1)/2)
		informative[x$nonfounders] = (ncols[x$nonfounders] < nall^2)
		for (i in 1:x$nInd) {
			if(ncols[i] == 1) next
			g = genolist[[i]]; kjonn = p[i, 'SEX']
			if(i %in% x$nonfounders && informative[far <- p[i, 'FID']])
				g = g[, g[1,] %in% genolist[[far]][1,] | g[1,] %in% genolist[[far]][2,], drop=F]
			if(i %in% x$nonfounders && informative[mor <- p[i, 'MID']])
				g = g[, g[2,] %in% genolist[[mor]][1,] | g[2,] %in% genolist[[mor]][2,], drop=F]
			barn = offs[[i]]
			for(b in barn[informative[barn]]) {
				g = g[, g[1,] %in% genolist[[b]][kjonn, ] | g[2,] %in% genolist[[b]][kjonn, ], drop=F]
			}
			genolist[[i]] = g
		}
		ncols_ny = unlist(lapply(genolist, ncol))
		if(any(ncols_ny == 0)) {attr(genolist, 'impossible') = TRUE; return(genolist)}
		if(sum(ncols_ny)==sum(ncols)) return(genolist)
	}
	genolist
}


.reduce_alleles = function(marker) {
	if(all(marker!=0)) return(marker)
	present = unique.default(as.numeric(marker))
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


#------------X-linked-------------------


.build_genolist_X <- function(x, marker, eliminate) {  #marker: marker matrix, dim = (nInd , 2). 
	n = attr(marker, 'nalleles')
	nseq = seq_len(n)
	.COMPLETE = { tmp1 = rep(nseq, each=n); tmp2 = rep.int(nseq, times=n);
		fath = c(tmp1, tmp2); moth = c(tmp2, tmp1)
		rbind(fath, moth, deparse.level=0)[, !duplicated.default(fath*1000 + moth), drop=F] #faster than unique(m, MARGIN=2)
	}
	SEX = x$pedigree[, 'SEX']; females = (1:x$nInd)[SEX==2]
	genolist = list()
	genolist[SEX==1 & marker[,1]==0] = list(nseq)
	genolist[SEX==1 & marker[,1]!=0] = marker[SEX==1 & marker[,1]!=0, 1]
	genolist[females] = lapply(females, function(i)  {
		g_i = marker[i, ]
		m = switch(sum(g_i == 0) + 1, 
		{ a=g_i[1]; b=g_i[2]; if(a==b) cbind(g_i, deparse.level=0) else cbind(g_i, c(b, a), deparse.level=0)},
		{ nz = g_i[g_i!=0]; rbind(c(nseq, rep.int(nz, n-1)), c(rep.int(nz, n), nseq[-nz]), deparse.level=0)},
		{ .COMPLETE })
		if ((i %in% x$founders) && (!x$orig.ids[i] %in% x$loop_breakers)) 
			m = m[, m[1,] <= m[2,], drop=FALSE]
		m
	})
	attr(genolist, 'impossible') = FALSE
	.eliminate_X(x, genolist, n, eliminate)
}



.eliminate_X = function(x, genolist, nall, repeats=0) {
	if(repeats==0 || attr(genolist, 'impossible')) return(genolist)
	SEX = x$pedigree[, 'SEX']; FID = x$pedigree[, 'FID']; MID = x$pedigree[, 'MID'] 
	males = (1:x$nInd)[SEX==1]; females = (1:x$nInd)[SEX==2]; fem_fou = .myintersect(females, x$founders); fem_nonfou = .myintersect(females, x$nonfounders)
	offs = lapply(1:x$nInd, function(i) offspring(x, i, original.id=FALSE))
	p = x$pedigree
	informative = logical(x$nInd)
	ncols_ny = unlist(lapply(genolist, length))/SEX #males are vectors, females matrices w/ 2 rows
	for(k in seq_len(repeats)) {
		ncols = ncols_ny
		informative[males] = (ncols[males] < nall)
		informative[fem_fou] = (ncols[fem_fou] < nall*(nall+1)/2)
		informative[fem_nonfou] = (ncols[fem_nonfou] < nall^2)
		for (i in males) {
			if(ncols[i] == 1) next
			g = genolist[[i]]
			if(i %in% x$nonfounders && informative[mor <- MID[i]]) 	g = g[ g %in% genolist[[mor]][1,] | g %in% genolist[[mor]][2,] ]
			barn = offs[[i]]
			for(b in barn[informative[barn] & SEX[barn]==2]) 	g = g[ g %in% genolist[[b]][1, ] ]
			genolist[[i]] = g
		}
		for (i in females) {
			if(ncols[i] == 1) next
			g = genolist[[i]]
			if(i %in% x$nonfounders && informative[far <- FID[i]])
				g = g[, g[1,] %in% genolist[[far]], drop=F]
			if(i %in% x$nonfounders && informative[mor <- MID[i]])
				g = g[, g[2,] %in% genolist[[mor]][1,] | g[2,] %in% genolist[[mor]][2,], drop=F]
			barn = offs[[i]]
			for(b in barn[informative[barn]]) {
				if(SEX[b]==1) g = g[, g[1,] %in% genolist[[b]] | g[2,] %in% genolist[[b]], drop=F]
				else g = g[, g[1,] %in% genolist[[b]][2, ] | g[2,] %in% genolist[[b]][2, ], drop=F]
			}
			genolist[[i]] = g
		}
		ncols_ny = unlist(lapply(genolist, length))/SEX
		if(any(ncols_ny == 0)) {attr(genolist, 'impossible') = TRUE; return(genolist)}
		if(sum(ncols_ny) == sum(ncols)) return(genolist)
	}
	genolist
}


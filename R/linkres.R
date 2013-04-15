print.linkres = function(x, ...) {
	rownames(x) = paste(rownames(x), ":", sep="")
   if(attr(x, 'analysis')=="mlink")
      rownames(x) = paste("theta=", rownames(x), sep="")
   df = as.data.frame(unclass(x))
   print(df, print.gap=2, ...)
}

summary.linkres = function(object, ...) {
	x = object
	switch( attr(x, "analysis"), 
	mlink = {
		mlods <- apply(x, 2, function(i) if (all(is.na(i))) NA else max(i, 0, na.rm=TRUE)); 
		MX <- max(mlods, na.rm=TRUE)
		cat("Max LOD score:", MX,"\n")
		cat("Achieved at marker(s):", names(mlods)[which(MX-mlods < 1e-4)],"\n")
	}, ilink = {
		lods=x['LOD',]; tmax=x['t_max',]; MX=max(lods, na.rm=TRUE)
		cat("Max LOD score:", MX, "\n")
		cat("Achieved at the following marker(s):\n")
		print(x[, which(MX - lods < 1e-4), drop=FALSE])
	})
}

as.data.frame.linkres = function(x, ..., sort=TRUE) {
   map = attr(x, 'map')
   map = map[map$MARKER %in% colnames(x), ,drop=FALSE]
   if(sort) map = map[order(map$CHR, map$POS), , drop=FALSE]
   if(attr(x,'analysis')=="mlink") 
      LOD = apply(x[, map$MARKER, drop=FALSE], 2, function(co) if(all(is.na(co))) NA else max(co, na.rm=TRUE))
   else {LOD = x['LOD',]; x = x[2,,drop=FALSE]}
	lods = cbind(map, LOD=LOD, t(x[, map$MARKER, drop=FALSE]))
   lods
}

lod.peaks = function(x, threshold, width=1) { # x et linkres objekt, eller data.frame med CHR, POS, LOD
   
   peak1chr = function(xchr, threshold, width) { # x must be sorted!
      rl = rle(xchr$LOD >= threshold)
      while(1) {
         short = rl$values & (rl$lengths < width)
         if(any(short)) {
            rl$values[short] <- FALSE
            rl = rle(inverse.rle(rl))
         } else break
      }
      
      xchr_nrow = nrow(xchr)
      if(!any(rl$values)) return(NULL)
      start_ind = c(0, cumsum(rl$lengths))[which(rl$values)]
      stop_ind = start_ind + rl$lengths[rl$values] + 1 # plus 1 to compensaate for endpoint[1]
      L = length(start_ind)
      if(start_ind[1]==0)       stop_ind[1] = stop_ind[1]-1
      if(stop_ind[L] > xchr_nrow) stop_ind[L] = xchr_nrow
      lapply(1:L, function(i) {p = xchr[start_ind[i]:stop_ind[i], ,drop=F]; rownames(p)=NULL; p})
   }
   
   
   df = as.data.frame(x)
   df = df[!is.na(df$LOD), ]
   chrs = unique.default(df$CHR)
   res = list()
   for(chr in chrs) {
      dfchr = df[df$CHR==chr, ]
      res = c(res, peak1chr(dfchr, threshold=threshold, width=width))
   }
   res
}    



.getMap = function(x, markers=seq_len(x$nMark), na.action=0, verbose=TRUE) {
	m=x$markerdata[markers]
	chrom=unlist(lapply(m, attr, 'chrom')); marker=unlist(lapply(m, attr, 'name')); pos=unlist(lapply(m, attr, 'pos'))
	map = data.frame(CHR=chrom, MARKER=marker, POS=pos, stringsAsFactors=FALSE)
	if(na.action > 0) {
		nas = is.na(marker) 	
		map$MARKER[nas] = paste("M", markers[nas], sep="")
	}
	if(na.action == 1) {
		nas2 = (is.na(chrom) | is.na(pos))
		if (all(nas2)) {
			if(verbose) cat("Warning: No map info given. Creating dummy map.\n")
			map$CHR = rep.int(1, x$nMark)
			map$POS = seq_len(x$nMark)
		}
	}
	# if(na.action ==2) 
		# nas2 = (is.na(chrom) | is.na(pos))
		# if(any(nas2)) {
			# if(verbose) cat("Warning: Deleting", sum(nas2), "markers with missing map coordinates.\n")
			# map = map[!nas2, , drop=FALSE]
		# }
	map
}

plot.linkres=function(x, chrom=NULL, ylim=NULL, ...) {
	analysis = attr(x, "analysis")
	map = attr(x, 'map')
	if(any(is.na(map$CHR))) stop("Incomplete or missing map.")
	
	map = map[map$MARKER %in% colnames(x), ,drop=FALSE]
	map = map[order(map$CHR, map$POS), , drop=FALSE]
	x = x[, match(map$MARKER, colnames(x)), drop=FALSE]
	
	if(!is.null(chrom)) {
		subindex = which(map$CHR %in% chrom)
		if(length(subindex)==0) stop("No markers on indicated chromosome(s).")
		subx = structure(x[, subindex, drop=FALSE], analysis=analysis, map=map[subindex, ,drop=F], class="linkres")
		return(plot(subx, chrom=NULL, ylim=ylim, ...))
	}
	
	switch(analysis, 
		mlink = {lds <- apply(x,2,max); if (is.null(ylim)) ylim <- c(-1.2, max(c(3, lds), na.rm=T) + 0.3)}, 
		ilink = {lds <- x['LOD',]; if (is.null(ylim)) ylim <- c(-0.5, max(c(3, lds), na.rm=T) + 0.3)},
	)

	nM = ncol(x)
	pos = map$POS
	chr_br = which(map$CHR[-1] != map$CHR[-nM])
	for(b in chr_br) pos[(b+1):nM] = pos[(b+1):nM] + map$POS[b] # NB: by now, map$POS != pos
	
	multichr = length(chr_br) > 0
	
	plot(pos, sapply(lds, max, ylim[1]), ylim=ylim, type="l", lwd=2, cex=.3, 
		xlab=ifelse(multichr, "Chromosome", paste("Position (cM) on chromosome", map$CHR[1])), 
		xaxt=ifelse(multichr, "n", par("xaxt")), 
		ylab="LOD score", ...)
	
	if (multichr) axis(1, at=c(0,pos[chr_br]), labels = map$CHR[c(chr_br, nM)], lwd.ticks=2)
}

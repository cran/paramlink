subset.linkdat <- function(x, subset=1:x$nInd, ..., markers=seq_len(x$nMark)) {
	#a few checks:
	for (i in subset) {
		offs=offspring(x,i)
		if (length(offs)==0) next
		if ( any( !(offs %in% subset) &  sapply(offs, function(b) any(offspring(x,b) %in% subset) ) ) )
			stop(paste("Individual", i,"has grandchildren in the subset, but no in-betweens are included." ))
	}
	xframe=as.data.frame(x, markers=markers)
	
	newfr <- subset(xframe, subset = xframe[,'ID'] %in% subset)
	newfr[!(newfr[,'FID'] %in% newfr[,'ID']), 'FID'] <- 0  #set FID=0 if father is not in subset
	newfr[!(newfr[,'MID'] %in% newfr[,'ID']), 'MID'] <- 0  #set MID=0 if mother is not in subset
 
	mis = ifelse(is.null(x$markerdata), 0, attr(x$markerdata, "missing"))
	linkdat(newfr, model=x$model, missing=mis)
}

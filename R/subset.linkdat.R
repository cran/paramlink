subset.linkdat <- function(x, subset=x$orig.ids, ..., markers=seq_len(x$nMark)) {
	x = removeMarkers(x, setdiff(seq_len(x$nMark), markers))
	xframe = as.matrix(x)
	
	newfr <- subset(xframe, subset = xframe[, 'ID'] %in% subset, ...)
	newfr[!(newfr[, 'FID'] %in% newfr[, 'ID']), 'FID'] <- 0  #set FID=0 if father is not in subset
	newfr[!(newfr[, 'MID'] %in% newfr[, 'ID']), 'MID'] <- 0  #set MID=0 if mother is not in subset
 
	restore_linkdat(newfr, attributes(xframe))
}

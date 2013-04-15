setMap = function(x, map, dat, pos=NULL, verbose=TRUE) {
	if(x$nMark==0) {x$map=NULL; return(x)}
	if(!is.null(pos)) {
		if(length(pos) != x$nMark) stop("Length of 'pos' argument does not match the number of markers.")
		x$map = data.frame(CHR=1, MARKER=paste("M", 1:x$nMark, sep=""), POS=pos, stringsAsFactors=FALSE)
		return(x)
	}
	
	stopifnot(!missing(map), !is.null(map), (is.data.frame(map) || is.character(map)))
	if(is.data.frame(map)) {
		if(ncol(map)>=3 && nrow(map)==x$nMark) {names(map)[1:3] = c('CHR', 'MARKER', 'POS'); x$map = map}
		else warning("Map not set: Something is wrong with the 'map' data frame.")
		return(x)
	}
	
	stopifnot(!missing(dat), !is.null(dat), is.character(dat))
	rawmap = read.table(map, as.is=TRUE, header=FALSE)
	if(!any(1:9 %in% strsplit(rawmap[1,1],"")[[1]])) rawmap = rawmap[-1,] #If no number occurs in first entry, first row is assumed to be header. 
	rawmap[[1]][rawmap[[1]]=="X"] = 23
	map1 = data.frame(CHR = as.numeric(rawmap[,1]), MARKER = as.character(rawmap[, 2]), POS = as.numeric(rawmap[, 3]), stringsAsFactors=FALSE) 
	rawdat = read.table(dat, as.is=TRUE, header=FALSE)
	dat = as.character(rawdat[rawdat[1]=="M", 2])  #names of all markers in map
		
	Mmatch = match(dat, map1$MARKER, nomatch=0)
	if(any(Mmatch == 0)) {
		del = dat[Mmatch==0]
		if(verbose) cat("Deleting the following marker(s), which are not found in the map file:\n", paste(del, collapse="\n"), "\n")
		x$markerdata[Mmatch==0] = NULL
	}

	map = map1[Mmatch, ]
	map = map[order(map$CHR, map$POS), ] 
		
	ord = match(dat, map$MARKER, nomatch=0)
	x$markerdata = x$markerdata[ord]
	x$nMark = length(x$markerdata)
	
	x$map = map
	x
}

.readDatfile = function(datfile, chrom) {
	dat = apply(read.table(datfile, as.is=T, comment.char="<", fill=T), 1, function(r) as.character(r[!is.na(r) & nzchar(r)]))
	nMark = as.numeric(dat[[1]])[1] - 1
	stopifnot(all(as.numeric(dat[[3]])==seq_len(nMark+1)), length(dat)==7+2*nMark+3)
	
	markernames = sapply(dat[6+(1:nMark)*2], '[', 4)
	pos = cumsum(as.numeric(dat[[length(dat)-1]]))
	equal = (pos[-1]==pos[-length(pos)]); k=0
	for(i in 2:length(pos)) 
		if (equal[i-1]) pos[i]=pos[i] + 0.0001*(k <- k+1) else k=0  #if consecutive entries are equal, add 0.0001's. 
		
	map = data.frame(CHR=chrom, MARKER=markernames, POS=pos, stringsAsFactors=F)
	
	freqlist = lapply(dat[7+(1:nMark)*2], function(r) as.numeric(r))
	names(freqlist) = markernames
	
	list(map=map, freq=freqlist)
}
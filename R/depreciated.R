# don't remember where this was used
.setMap = function(x, map, dat, pos = NULL, verbose = TRUE) {
    if (x$nMark == 0) {
        x$map = NULL
        return(x)
    }
    if (!is.null(pos)) {
        if (length(pos) != x$nMark) 
            stop("Length of 'pos' argument does not match the number of markers.")
        x$map = data.frame(CHR = 1, MARKER = paste("M", 1:x$nMark, sep = ""), POS = pos, stringsAsFactors = FALSE)
        return(x)
    }
    
    stopifnot(!missing(map), !is.null(map), (is.data.frame(map) || is.character(map)))
    if (is.data.frame(map)) {
        if (ncol(map) >= 3 && nrow(map) == x$nMark) {
            names(map)[1:3] = c("CHR", "MARKER", "POS")
            x$map = map
        } else warning("Map not set: Something is wrong with the 'map' data frame.")
        return(x)
    }
    
    stopifnot(!missing(dat), !is.null(dat), is.character(dat))
    rawmap = read.table(map, as.is = TRUE, header = FALSE)
    if (!any(1:9 %in% strsplit(rawmap[1, 1], "")[[1]])) 
        rawmap = rawmap[-1, ]  #If no number occurs in first entry, first row is assumed to be header. 
    rawmap[[1]][rawmap[[1]] == "X"] = 23
    map1 = data.frame(CHR = as.numeric(rawmap[, 1]), MARKER = as.character(rawmap[, 2]), POS = as.numeric(rawmap[, 
        3]), stringsAsFactors = FALSE)
    rawdat = read.table(dat, as.is = TRUE, header = FALSE)
    dat = as.character(rawdat[rawdat[1] == "M", 2])  #names of all markers in map
    
    Mmatch = match(dat, map1$MARKER, nomatch = 0)
    if (any(Mmatch == 0)) {
        del = dat[Mmatch == 0]
        if (verbose) 
            cat("Deleting the following marker(s), which are not found in the map file:\n", 
                paste(del, collapse = "\n"), "\n")
        x$markerdata[Mmatch == 0] = NULL
    }
    
    map = map1[Mmatch, ]
    map = map[order(map$CHR, map$POS), ]
    
    ord = match(dat, map$MARKER, nomatch = 0)
    x$markerdata = x$markerdata[ord]
    x$nMark = length(x$markerdata)
    
    x$map = map
    x
}


# Not used - remove?
.SNPfreq = function(markernames, Bfreq, allele1 = "1", allele2 = "2", file = NULL) {
    # Create and write freq file in extended Merlin format.  Bfreq numerical vector with same
    # length as markernames (contains freqs for allele 2)
    
    stopifnot(length(markernames) == length(Bfreq))
    n = length(markernames)
    col1 = rep(c("M", "A", "A"), n)
    col2 = as.character(rbind(markernames, allele1, allele2))
    col3 = as.character(rbind("", 1 - Bfreq, Bfreq))
    res = cbind(col1, col2, col3)
    if (!is.null(file)) 
        write.table(res, file = file, col.names = F, row.names = F, quote = F)
    invisible(res)
}

.my.grid = fast.grid
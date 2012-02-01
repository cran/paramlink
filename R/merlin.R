merlin = function(x, markers=seq_len(x$nMark), model=TRUE, theta=NULL, options="", verbose=FALSE, generate.files=TRUE, cleanup=generate.files, outputfile="") {
	
	clean = function(cleanup, verbose, files) if (cleanup) {unlink(files); if(verbose) cat("Files successfully removed\n")}
	
	if (x$nMark == 0) stop("No markers exist for this linkdat object.")
	x = removeMarkers(x, seq_len(x$nMark)[-markers])
	map = paramlink:::.getMap(x, na.action=1, verbose=F)
	mNames = map$MARKER
	
	extensions = c("ped", "dat", "map", "freq", if(model) "model")
	
	if (generate.files) {
		files = write.linkdat(x, prefix="merlin", what=extensions, merlin=TRUE)
		if (verbose) {cat("Files successfully generated:\n");  print(files)}
	}	
	
	options = paste(options, "--markerNames --quiet")
	
	if (nonz_theta <- any(theta > 0)) {
		if (length(markers) > 1) {
			clean(cleanup, verbose, files)
			stop("Nonzero 'theta' values are possible only with a single marker.")
		}
		pos = as.numeric(map[1, 3]) - 50 * log(1 - 2 * theta) #Haldane's map: Converting rec.fractions to cM positions.
		options = paste(options, " --positions:", paste(pos, collapse=","), sep="")
	}
	
	command = paste("merlin -p merlin.ped -d merlin.dat -m merlin.map -f merlin.freq ", 
					if(model) "--model merlin.model ", options, sep="")
					
	if (verbose) cat("\nExecuting the following command:\n", command, "\n\n", sep="")

	merlinout = system(command, intern=T)
	clean(cleanup, verbose, files)
	if (nzchar(outputfile))		write(merlinout, outputfile)
	if (any(substr(merlinout,1,11) == "FATAL ERROR")) stop(paste(merlinout, collapse="\n"))
	if(verbose) {cat("Merlin run completed\n"); print(merlinout)}
	if (!is.na(skipped <- which(substr(merlinout,3,9) == "SKIPPED")[1])) stop(paste(merlinout[c(skipped-1, skipped)], collapse="\n"))
	
	if(!model) return(merlinout)

	##Extract LOD scores
	nchars = nchar(merlinout)
	chromStarts = which(substr(merlinout, 1, 20) == "Analysing Chromosome")
	chroms = as.numeric(substr(merlinout[chromStarts], 22, nchars[chromStarts]))
	if(length(chromStarts)==0) {
		chromStarts = match("Parametric Analysis", substr(merlinout, 1, 19))
		chroms = map$CHR[1]
	}
	
	lods.list = list()
	for(i in seq_along(chromStarts)) {
		outp = merlinout[-(1:chromStarts[i])]
		LODstart = 1 + match("POSITION", substr(outp,8,15))
		LODstop = LODstart -1 + match("", outp[-(1:LODstart)])
		
		lods.list = c(lods.list, lapply(outp[LODstart:LODstop], function(char) 
			c(chroms[i], scan(textConnection(char), what = "", quiet=TRUE)[1:2])))
	}
	lods.df = do.call(rbind, lods.list)
	
	mlods = lods.df[, 3]
	mlods[mlods == "-INFINITY"] = -Inf

	if(nonz_theta) {
		mlodsdim = c(length(theta), 1)
		dimnam = list(theta, mNames)
	}
	else {
		if(lods.df[1,2] %in% mNames) 
			markernames = lods.df[, 2]
		else {
			markernames = paste(lods.df[, 1], lods.df[, 2], sep="_") #if first "POSITION" entry is not a marker name, create new map with markernames formed as "chr_pos".
			map = data.frame(CHR=as.numeric(lods.df[, 1]), MARKER = markernames, POS = as.numeric(lods.df[, 2]))
		}
		mlodsdim = c(1, length(mlods))
		dimnam = list(0, markernames)
	}

	res = structure(as.numeric(mlods), dim=mlodsdim, dimnames=dimnam, analysis="mlink", map=map, class="linkres")
	res
}

# .superlink = function(x, theta=0) {
	# write.table(as.data.frame(x, famid=TRUE, missing=0), file="superlink.preR", col.names=F, row.names=F, quote=F)
	# gap:::makeped("superlink.preR", "superlink.pedR", 1)
	
	# cat("2 0 0 5 0 \n0 0.0 0.0 0 \n1 2 \n1 2\n", 
	# paste(format(c(1-x$model$dfreq, x$model$dfreq), scientific=F, decimal.mark="."), collapse=" "),
	# "\n1\n", paste(format(x$model$penetrances, scientific=F, decimal.mark="."), collapse=" "),
	# "\n3 ", x$model$nallel, "\n", paste(format(as.numeric(x$model$afreq), scientific=F, decimal.mark="."), collapse=" "),
	# "\n0 0\n", paste(format(theta, scientific=F, decimal.mark="."), collapse=" "),
	# "\n1 .51 0.5\n", sep="", file="superlink.datR")
	
	# super = system("superlink superlink.datR superlink.pedR", intern=T)
	# file.remove(c("superlink.preR","superlink.datR","superlink.pedR"))
	# super
# }

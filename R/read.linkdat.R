read.linkdat <-
function(file, header = FALSE, missing=0, model=NULL, ...) {
	dat <- read.table(file, header=header, as.is=TRUE, ...)
	# if ( (header && toupper(colnames(dat)[1])=="FAMID") || length(unique(dat[,1]))==1 ) {
		# cat("Removing first column (family id)\n")
		# dat=dat[,-1] 
		# if (!is.na(simcol)) simcol=simcol-1
	# }
	if (!header) {
		colnames(dat)[1:5] <- c('ID','FID','MID','SEX','AFF')
		# if (!is.na(simcol)) {
			# stopifnot(is.numeric(simcol), length(simcol)==1, simcol>5) 
			# colnames(dat)[simcol] <- 'SIM'
		# }
	}
	linkdat(dat, model=model, missing=missing)
}


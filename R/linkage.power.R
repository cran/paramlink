linkage.power <- function(x, N=100, all=FALSE, loop_breakers=NULL, threshold=NULL, seed=NULL) {
	if (is.null(x$model)) stop("No model set.")
	if (any(!is.numeric(N), length(N)>1, N%%1 != 0)) stop("N must be a positive integer.")
	
	if (all || length(x$available)==0) x = setAvailable(x, x$orig.ids)
	 
	sims = SNPsim(x, N=N, partialmarker=NULL, available=x$available, loop_breakers=loop_breakers, unique=FALSE, seed=seed)
	cat(N, "markers simulated\n")
	lds = structure( lod(sims, theta=0, loop_breakers=loop_breakers, verbose=FALSE), analysis="power", class="linkres")
	summary(lds, threshold=threshold)
	invisible(lds)
}


power.varyPar <- function(x, N=100, varyPar, values, all=FALSE, loop_breakers=NULL, seed=NULL) {
	if (is.null(x$model)) stop("No model set.")
	if (any(!is.numeric(N), length(N)>1, N%%1 != 0)) stop("N must be a positive integer.")
	
	if (all) x = setAvailable(x, x$orig.ids)
	
	models = switch(x$model$chr,
		AUTOSOMAL = {switch(varyPar,
			f0 = lapply(values, function(k) {mod=unclass(x$model); mod$penetrances[1]=k; mod}),
			f1 = lapply(values, function(k) {mod=unclass(x$model); mod$penetrances[2]=k; mod}),
			f2 = lapply(values, function(k) {mod=unclass(x$model); mod$penetrances[3]=k; mod}),
			dfreq = lapply(values, function(k) {mod=unclass(x$model); mod$dfreq=k; mod}),
			afreq = lapply(values, function(k) {mod=unclass(x$model); mod$afreq=c(k, 1-k); mod}),
			stop("Argument 'varyPar' must be one of 'f0', 'f1', 'f2', 'dfreq', 'afreq'.") ) },
		X = {switch(varyPar,
			f0_m = lapply(values, function(k) {mod=unclass(x$model); mod$penetrances$male[1]=k; mod}),
			f1_m = lapply(values, function(k) {mod=unclass(x$model); mod$penetrances$male[2]=k; mod}),
			f0_f = lapply(values, function(k) {mod=unclass(x$model); mod$penetrances$female[1]=k; mod}),
			f1_f = lapply(values, function(k) {mod=unclass(x$model); mod$penetrances$female[2]=k; mod}),
			f2_f = lapply(values, function(k) {mod=unclass(x$model); mod$penetrances$female[3]=k; mod}),
			dfreq = lapply(values, function(k) {mod=unclass(x$model); mod$dfreq=k; mod}),
			afreq = lapply(values, function(k) {mod=unclass(x$model); mod$afreq=c(k, 1-k); mod}),
			stop("Argument 'varyPar' must be one of 'f0_m', 'f1_m', 'f0_f', 'f1_f', 'f2_f', 'dfreq', 'afreq'.") ) } )

	sims <- SNPsim(x, N=N, loop_breakers=loop_breakers, seed=seed, unique=TRUE) 
	cat(N, "markers simulated;", ncol(sims$markerdata)/2, "unique.\n")

	lods = sapply(models, function(mod) max(lod(setModel(sims, model=mod), theta=0, loop_breakers=loop_breakers, verbose=FALSE)))
	plot(values, lods, type='l', xlab=varyPar, ylab='Maximum LOD')
	data.frame(value=values, LOD=lods)
}

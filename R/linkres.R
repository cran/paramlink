print.linkres <- function(x, ...) {
	print(as.data.frame(unclass(x)), ...)
}

summary.linkres <- function(object, threshold=NULL, ...) {
	x <- object
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
	}, power = {	
		cat("Highest singlepoint LOD score for simulated markers:", max(x),"\n")
		cat("Markers obtaining maximum score:", mean((max(x)-x) < 1e-4)*100, "%\n")
		cat("ELOD:", mean(x), "\n")
		if (!is.null(threshold)) {
			cat("\nThreshold:", threshold, "\n")
			cat("Markers with score above threshold:", mean(x>=threshold)*100, "%\n")
		}
	})
}


plot.linkres=function(x, ylim=NULL, ...) {
	switch(attr(x, "analysis"), 
		mlink = {lds <- apply(x,2,max); if (is.null(ylim)) ylim <- c(-1.2, max(c(3, lds), na.rm=T) + 0.3) }, 
		ilink = {lds <- x['LOD',]; if (is.null(ylim)) ylim <- c(-0.5, max(c(3, lds), na.rm=T) + 0.3)}, 
		power = {stop("Not implemented")}
	)
	plot(sapply(lds, max, ylim[1]), ylim=ylim, type="l", lwd=2, cex=.3, xlab="Markers", ylab="LOD score", ...)
}

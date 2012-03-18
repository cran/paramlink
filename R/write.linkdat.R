write.linkdat=function(x, file="", col.names=FALSE, famid=TRUE) { 
	df=as.data.frame(x)
	#if (use.oldid) {df['ID'] <- x$oldid; df['OLDID'] <- NULL}
	if (famid) df=cbind(FAMID=1, df)
	write.table(df, file=file ,sep="\t", quote=FALSE, col.names=col.names, row.names=FALSE)
}

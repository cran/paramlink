write.linkdat=function(x, file="", col.names=FALSE) {
   df = as.data.frame(x, famid=T)
   write.table(df, file=file, sep="\t", quote=FALSE, col.names=col.names, row.names=FALSE)
}

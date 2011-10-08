plot.linkdat <- function(x, marker=NULL, alleles=NULL, sep="", missing="-", id.labels=x$orig.ids, title.name=x$famid, margins=c(4.1, 1, 4.1, 1), ...) {

  ped <- x$pedigree;  AFF <- ped[,'AFF']; ID <- ped[,'ID']; SEX <- ped[,'SEX']
  cols <- ifelse(AFF==0,8,1)  #unknown affection tegnes grå
  AFF[AFF==0] <- 2 #må gjøre unknowns til affected for å unngå bug

  strid <- paste("", id.labels, sep="\n")

  if (!is.null(x$sim)) strid <- paste(strid, ifelse(x$sim==2, "*", ""), sep="")

  if (!is.null(marker)) {
	if (is.null(m <- x$markerdata) || marker > x$nMark) stop("Indicated marker does not exist")
	if (is.null(alleles)) alleles=attr(m, "alleles")[marker]
	else if(!is.list(alleles)) alleles=list(alleles)
	chrom = ifelse(is.null(x$model), 'AUTOSOMAL', x$model$chrom)
	geno = .prettyMarkers(m[, c(2*marker-1, 2*marker)], alleles=alleles, sep=sep, missing=missing, singleCol=TRUE, chrom=chrom, sex=SEX)
	strid <- paste("", strid, geno, sep="\n")
  }
  kped<-pedigree(id=ID, dadid=ped[,'FID'], momid=ped[,'MID'], sex=SEX, affected=AFF)
  par(xpd=T)
  plot.pedigree(kped, id=strid, col=cols, mar=margins, ...)
  title(paste("Family", title.name))
}



.haploplot <- function(x, haplotypes, id=FALSE, cex=1, margins=c(4.1,1,4.1,1), ...) {
  stopifnot(class(x)=="linkdat", is.matrix(haplotypes), (ncol(haplotypes)%%2)==0)
  
  ped <- x$pedigree;  AFF <- ped[,'AFF']; ID <- ped[,'ID']
  
  cols <- ifelse(AFF==0,8,1)  #unknown affection tegnes grå
  AFF[AFF==0] <- 2 #må gjøre unknowns til affected for å unngå bug

  N=ncol(haplotypes)/2
  all1 = haplotypes[, 2*seq_len(N)-1];   all2 = haplotypes[, 2*seq_len(N)] 
  geno = matrix(paste(all1, all2, sep="|"), ncol=N)
  geno = apply(geno, 1, paste, collapse="\n")
  
  if (id) strid <- paste("", ID, sep="\n\n") else strid <- rep("\n", x$nInd)
  strid <- paste("", strid, geno, sep="\n")
  
  if(is.null(mar)) mar=c(4.1+.5*N,2,4.1,2)
  kped<-pedigree(id=ID, dadid=ped[,'FID'], momid=ped[,'MID'], sex=ped[,'SEX'], affected=AFF)
  par(xpd=T)
  plot.pedigree(kped, id=strid, col=cols, cex=cex, mar=margins, ...)
}

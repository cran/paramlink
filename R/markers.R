setMarkers <- function(x, m, missing=0) {
	n = ncol(m)
	if (n%%2 != 0) stop("Uneven number of marker columns")
	m=as.matrix(m)
	alleles = sort(setdiff(unique(as.vector(m)), missing))
	if (length(alleles)==0) alleles <- c("1","2") 
	else if (length(alleles)==1) {
		if (alleles %in% c("1","2")) alleles <- c("1","2") 
		else if (alleles %in% c("a","b")) alleles <- c("a","b")
		else if (alleles %in% c("A","B")) alleles <- c("A","B")
		else if (alleles %in% c("x","y")) alleles <- c("x","y")
		else if (alleles %in% c("X","Y")) alleles <- c("X","Y")
	}
	mnum = match(m,alleles)
	dim(mnum) = dim(m)
	mnum[is.na(mnum)] = 0
	class(mnum) = "markerdata"
	attr(mnum,"alleles") <- alleles
	attr(mnum,"missing") <- missing
	x$markerdata = mnum; x$nMark = n/2
	x
}

getMarkers <- function(x, markers) {
	m = x[['markerdata']]
	.prettyMarkers(m[, sort.int(c(2*markers-1, 2*markers))], alleles=attr(m, "alleles"), missing=attr(m, "missing"), singleCol=FALSE)
}

.prettyMarkers = function(m, alleles=NULL, sep="", missing=NULL, mNames=NULL, singleCol=FALSE, chrom, sex) {
	if (is.null(m) || length(m)==0) return(m)
	n = ncol(m)/2
	if (is.null(alleles)) alleles=attr(m,"alleles")
	if (is.null(missing)) missing=attr(m,"missing")
	if (is.null(mNames)) mNames = paste('M', seq_len(n), sep="")

	orig.markers = c(missing, alleles)[m+1]
	dim(orig.markers) = dim(m)

	if (singleCol) {
		al1 = orig.markers[, 2*seq_len(n)-1, drop=FALSE]
		al2 = orig.markers[, 2*seq_len(n), drop=FALSE]
		if (missing(chrom) || is.null(chrom)) stop("Missing argument: 'chrom'.")
		switch(chrom, AUTOSOMAL = { 
			m.matrix = matrix(paste(al1, al2, sep=sep), ncol=n )
		}, X = {
			if (missing(sex)) stop("Missing argument: 'sex'.")
			males = which(sex==1); females = which(sex==2)
			if (!all(hh <- al1[males,] == al2[males,])) warning("Some male individuals are heterozygous.")
			m.matr.male = al1[males, , drop=FALSE]
			m.matr.female = matrix(paste(al1, al2, sep=sep), ncol=n)[females, , drop=FALSE]
			m.matrix = rbind(m.matr.male, m.matr.female)[order(c(males, females)), , drop=FALSE]
		})
		colnames(m.matrix) = mNames
		return(m.matrix)
	}
	else {
		colnames(orig.markers) = paste(rep(mNames, each=2), 1:2, sep="_")
		return(orig.markers)
	}
}

modifyMarker <- function(x, id, alleles, marker=1) {
	m = x[['markerdata']]
	if (!is.numeric(marker) || length(marker) != 1) stop("Argument 'marker' must be a single integer.")
	if (is.null(m) || marker > x[['nMark']]) stop("Indicated marker does not exist.")
	if (length(alleles)>2) stop("Argument 'alleles' must have length at most 2.")
	mis = attr(m, "missing")
	#if (!all(alleles %in% c(mis, attr(m, "alleles")))) stop("Valid alleles are", c(mis, attr(m, "alleles")))
	pm <- .prettyMarkers(m)
	for (i in id) pm[i, c(2*marker-1, 2*marker)] <- alleles
	setMarkers(x, pm, missing=mis)
}

.diallel2geno <- function(marker) { #marker a numerical nInd * 2 matrix   
	#Coding genotypes as single integer: 00 -> 0, 11 -> 1, 22 -> 2, 12/21 -> 3, 01/10 -> 4, 02/20 -> 5.
	#Each pair of alleles is seen as an integer written in base 3, and this integer is permuted to fit with the above code.
	c(0,4,5,4,1,3,5,3,2)[colSums(c(3,1) * t(marker)) + 1]
}

.geno2diallel <- function(codedgenos) { #input: matrix of single-numerical genotypes. 
#Ouput: Matrix with twice the number of columns, decoded as 1 -> 1 1, 2 -> 1 2, 3 -> 2 2
decode = sapply(t(codedgenos+1),function(i) switch(i, c(0,0), c(1,1), c(2,2), c(1,2)))
dim(decode) = c(2*ncol(codedgenos), nrow(codedgenos))
t(decode)
}

mendelianCheck <- function(x) {
  getAlleles = function(m, id, marker) as.vector(m[id, c(2*marker-1,2*marker)])  

  autosCheck = function(al.fa, al.mo, al.of) {
	fa0 = (0 %in% al.fa); mo0 = (0 %in% al.mo); of0 = (al.of==0)
	ff = fa0 | of0 | (al.of %in% al.fa)
	mm = mo0 | of0 | (al.of %in% al.mo)
	(ff[1] && mm[2]) || (ff[2] && mm[1])
  }
  
  maleXHomoz = function(al.of) 
	al.of[1]==al.of[2]
	
  maleXCheck = function(al.mo, al.of)
	(al.of==0) || (0 %in% al.mo) || (al.of %in% al.mo)
	
  p = x$pedigree
  m = x$markerdata
  nucs = x$subnucs
  chr = ifelse(is.null(x$model), "AUTOSOMAL", x$model$chrom)
  
  error = numeric(0)
  switch(chr, 
  AUTOSOMAL = {
	for (i in seq_along(x$nMark)) for (nuc in nucs) {
		al.fa = getAlleles(m, nuc[2], i)
		al.mo = getAlleles(m, nuc[3], i)
		for (of in nuc[-c(1:3)])
			if (!autosCheck(al.fa, al.mo, getAlleles(m, of, i))) {
			  cat("Individual ", of, ", marker ",i, ": Alleles not compatible with parents.\n", sep="")
			  error = union(error, i)
			}
	}
  }, 
  X = {
    sex = p[, 'SEX']
	for (i in seq_along(x$nMark)) for (nuc in nucs) {
		al.fa = getAlleles(m, nuc[2], i)
		al.mo = getAlleles(m, nuc[3], i)
		for (of in nuc[-c(1:3)]) {
			al.of = getAlleles(m, of, i)
			switch(sex[of], {
				if (!maleXHomoz(al.of)) {
					cat("Individual ", of, ", marker ",i, ": Male heterozygosity not compatible with X-linked model.\n", sep="")
					error = union(error, i)
				}
				if (!maleXCheck(al.mo, al.of)) {
					cat("Individual ", of, ", marker ",i, ": Allele not compatible with mother.\n", sep="")
					error = union(error, i)
				}
			}, {
			  if (!autosCheck(al.fa, al.mo, getAlleles(m, of, i))) {
			    cat("Individual ", of, ", marker ",i, ": Alleles not compatible with parents.\n", sep="")
			    error = union(error, i)
			  }
			})
		}
	}
  })
  error
}
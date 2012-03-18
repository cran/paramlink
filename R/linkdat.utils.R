offspring <- function(x,i) {
	p=x[['pedigree']]
	sex=p[i,'SEX']
	which(as.numeric(p[, sex+1])==i)
}	


spouses <- function(x, i) {
	p=x$pedigree
	sex=p[i, 'SEX']
	offs=offspring(x, i)
	unique(as.numeric(p[offs, 4-sex]))  #sex=1 -> column 3; sex=2 -> column 2.
}	

descendents=function(x,i) {
	nextgen <- desc <- offspring(x,i)
	while(1) {
    	nextgen <- unique(unlist( sapply(nextgen,offspring,x=x), use.names=F))
		if (length(nextgen)==0) 
			break
		desc <- union(desc,nextgen)
	}
	desc
}

ped.dist.matrix <- function(x) { 
	comb2 = function(n) {
		if (n<2) return(NULL)
		v1 = rep.int(seq_len(n-1), (n-1):1)
		v2 = NULL
		for (i in 2:n) v2 = c(v2, i:n)
		cbind(v1,v2, deparse.level=0)
	} 

	d = matrix(0, ncol=x$nInd, nrow=x$nInd)
	p=x$pedigree
	for (f in 1:x$nInd) {
		d[ f, offspring(x,f) ] <- 1
		d[ f, p[f, c('FID', 'MID')] ] <- 1
	}
	k=1
	while (prod(d[col(d)!=row(d)])==0) { 
		k=k+1
		for (r in 1:nrow(cc <- comb2(x$nInd))) {
			a=cc[r,1]; b=cc[r,2]
			if (d[a,b] > 0) next
			if (any( d[a,]>0 & d[,b]>0 & d[a,]+d[,b] == k))
			d[a,b] <- d[b,a] <- k
		}
	}
	d
}

getpaths = function(obj, i, j) {
	if (class(obj)=="linkdat") d=ped.dist.matrix(obj)
	else if (is.matrix(obj)) d=obj
	else stop("'obj' must be either a matrix or a 'linkdat' object.")
	if (i==j) return(i)
	D=d[i,j]
	res=matrix(i)
	for (k in 1:D) {
		tmp=list()
		for (r in 1:nrow(res)) {
			laststep = res[r, k]
			newsteps=which( d[laststep, ] == 1 & d[, j] == D-k)
			tmp[[r]] = cbind( matrix(res[r, 1:k], ncol=k, nrow=length(newsteps),byrow=TRUE), newsteps, deparse.level=0)
		}
		res = do.call(rbind, tmp)
	}
	res
}

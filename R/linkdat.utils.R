.internalID = function(x, orig.ids) {
	internal_ids = match(orig.ids, x$orig.ids)
	if(any(is.na(internal_ids))) stop(paste("Indicated ID(s) not among original ID labels:", paste(orig.ids[is.na(internal_ids)], collapse=",")))
	internal_ids
}
.getSex = function(x, ids) {as.vector(x$pedigree[match(ids, x$orig.ids), 'SEX'])}
.comb2 = function(n) {
		if (n<2) return(matrix(nrow=0, ncol=2))
		v1 = rep.int(seq_len(n-1), (n-1):1)
		v2 = NULL
		for (i in 2:n) v2 = c(v2, i:n)
		cbind(v1,v2, deparse.level=0)
	} 

offspring <- function(x, id, original.id=TRUE) {
	internal_id = ifelse(original.id, .internalID(x, id), id)
	p = x$pedigree
	offs_rows = p[, 1 + p[internal_id, 'SEX'] ] == internal_id
	if(original.id) return(x$orig.ids[offs_rows]) else return(which(offs_rows))
}	


spouses <- function(x, id, original.id=TRUE) { #returns a vector containing all individuals sharing offspring with 'i'.
	internal_id = ifelse(original.id, .internalID(x, id), id)
	p = x$pedigree
	offs_rows = p[, 1 + p[internal_id, 'SEX'] ] == internal_id
	spou = unique(p[offs_rows, 4 - p[internal_id, 'SEX']])  #sex=1 -> column 3; sex=2 -> column 2.
	if (original.id) return(x$orig.ids[spou]) else return(spou)
}	


descendants=function(x, id, original.id=TRUE) { 
	internal_id = ifelse(original.id, .internalID(x, id), id)
	nextgen <- desc <- offspring(x, internal_id, original.id=FALSE)
	while(TRUE) {
    	nextgen <- unlist( lapply(nextgen, offspring, x=x, original.id=FALSE))
		if (length(nextgen)==0) break
		desc <- c(desc, nextgen)
	}
	desc = unique(sort.default(desc))
	if (original.id) return(x$orig.ids[desc]) else return(desc)
}


#pedDistMatrix(x) computes a symmetrix integer matrix (d_ij), where d_ij is the pedigree distance (i.e., the length of the shortest pedigree path) between i and j. For example, the following relations have the indicated pedigree distances: parent/offspring	1, siblings	2, parents	2, uncle/niece	3, first cousins 4.
.pedDistMatrix <- function(x) { 
	d = matrix(0, ncol=x$nInd, nrow=x$nInd)
	for (f in 1:x$nInd) {
		d[ f, offspring(x, f, original.id=FALSE) ] <- 1
		d[ f, x$pedigree[f, c('FID', 'MID')] ] <- 1
	}
	k=1
	while (prod(d[col(d)!=row(d)])==0) { 
		k=k+1
		for (r in 1:nrow(cc <- .comb2(x$nInd))) {
			a=cc[r,1]; b=cc[r,2]
			if (d[a,b] > 0) next
			if (any( d[a,]>0 & d[,b]>0 & d[a,]+d[,b] == k))
			d[a,b] <- d[b,a] <- k
		}
	}
	d
}

.getpaths = function(obj, i, j) {
	if (class(obj)=="linkdat") d = .pedDistMatrix(obj)
	else if (is.matrix(obj)) d = obj
	else stop("'obj' must be either a matrix or a 'linkdat' object.")
	if (i==j) return(i)
	D = d[i,j]
	res = matrix(i)
	for (k in 1:D) {
		tmp = list()
		for (r in 1:nrow(res)) {
			laststep = res[r, k]
			newsteps = which( d[laststep, ] == 1 & d[, j] == D-k)
			tmp[[r]] = cbind( matrix(res[r, 1:k], ncol=k, nrow=length(newsteps),byrow=TRUE), newsteps, deparse.level=0)
		}
		res = do.call(rbind, tmp)
	}
	res
}

pedigreeLoops = function(x) {
	dls = .descentPaths(x, 1:x$nInd, original.ids=FALSE)
	loops = list()
	for (id in 1:x$nInd) {
		if (length(dl <- dls[[id]])==1) next
		pairs = .comb2(length(dl))
		for (p in 1:nrow(pairs)) {
			p1 = dl[[pairs[p,1]]]; p2 = dl[[pairs[p,2]]]
			if (p1[2]==p2[2]) next
			inters = p1[match(p2,p1,0L)][-1]  #intersecting path members, excluding id
			if (length(inters)==0) next
			else {
				top=x$orig.ids[p1[1]]; bottom=x$orig.ids[inters[1]]
				pathA = p1[ seq_len(which(p1==inters[1]) - 2) + 1 ] #without top and bottom. Seq_len to deal with the 1:0 problem.
				pathB = p2[ seq_len(which(p2==inters[1]) - 2) + 1 ]
				loops = c(loops, list(list(top=top, bottom=bottom, pathA = x$orig.ids[pathA], pathB=x$orig.ids[pathB])))
			}
		}
	}
	unique(loops)
}

.descentPaths = function(x, ids, original.ids = TRUE)  {
	if (original.ids) ids = .internalID(x, ids)
	offs = lapply(1:x$nInd, offspring, x=x, original.id=FALSE)
	lapply(ids, function(id) {
		res=list(id)
		while(TRUE) {
			newoffs = offs[sapply(res, function(path) path[length(path)])]
			if (length(unlist(newoffs))==0) break

			nextstep = lapply(1:length(res), function(r) if   (length(newoffs[[r]])==0) res[r]   else lapply(newoffs[[r]], function(kid) c(res[[r]], kid) ))
			res = unlist(nextstep, recursive=FALSE)
		}
		if (original.ids) lapply(res, function(internal_vec) x$orig.ids[internal_vec]) else res
	})
}


.siblings = function(x, id, inclusive=FALSE, original.id=TRUE) {
	if(original.id) id = .internalID(x, id)
	if(id %in% x$founders) sibs = id
	else {
		p = x$pedigree
		sibs = which((p[, 'FID'] == p[id, 'FID']) & (p[, 'MID'] == p[id, 'MID']))
	}
	if(!inclusive) sibs = setdiff(sibs, id)
	if(original.id) sibs = x$orig.ids[sibs]
	sibs
}

uninformative = function(x, return.reduced.ped=FALSE) {
	p = relabel(x$pedigree, x$orig.ids)
	remove_IDs = function(p, ids) p[!p[,'ID'] %in% ids, ]
	
	uninformative = numeric(0)
	while(TRUE) {
		ID = p[, 'ID']
		leaves = setdiff(ID, p[, c('FID','MID')])
		nonaff = p[ p[,'AFF']<2 , 'ID' ]; healthy = p[ p[,'AFF']==1 , 'ID' ]; unknown = p[ p[,'AFF']==0 , 'ID' ]
		uninf = intersect(leaves, unknown) #these are always uninformative

		healthy_leaves = intersect(leaves, healthy)
		nonaff_parents = (p[match(healthy_leaves, ID), 'FID'] %in% nonaff) & (p[match(healthy_leaves, ID), 'MID'] %in% nonaff)
		maybe = healthy_leaves[nonaff_parents] #these are healthy leaves with healthy or unknown parents. Will be marked as uninf iff all sibs are healthy
		if(length(maybe) > 0) {
			healthy_sibs = sapply(maybe, function(id) {intid = match(id, ID); sibs = p[ (p[, 'FID'] == p[intid, 'FID']) & (p[, 'MID'] == p[intid, 'MID']), 'ID']; all(sibs %in% healthy)})
			uninf = c(uninf, maybe[healthy_sibs])
		}
		if(length(uninf)==0) break
		
		p = remove_IDs(p, uninf)
		fou = p[ p[, 'FID']==0 & p[, 'MID']==0, 'ID']
		leftover_spouses = setdiff(fou, p[, c('FID','MID')])
		p = remove_IDs(p, leftover_spouses)
		uninformative = c(uninformative, uninf, leftover_spouses)
	}
	if(return.reduced.ped) {
		cat("Removing individuals:", paste(sort(uninformative), collapse=", "), "\n") 
		return(linkdat(p, model=x$model, verbose=FALSE))
	}
	else return(sort(uninformative))
}

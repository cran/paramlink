addMarker = function(x, m, alleles=NULL, afreq=NULL, missing=0) {
	if(is.matrix(m) || is.data.frame(m)) stopifnot(nrow(m)==x$nInd, ncol(m)==2)
	if(class(m)=="marker") 
		m = list(m)
	if(is.list(m) && all(unlist(lapply(m, class))=="marker")) 
		return(setMarkers(x, structure(c(x$markerdata, m), class="markerdata")))
	if(!is.list(m) && length(m)==1)	
		m = matrix(m, ncol=2, nrow=x$nInd)  #gives a nice way of setting an empty or everywhere-homozygous marker, e.g.: x=addMarker(x,0)
	mm = .createMarkerObject(m, alleles, afreq=afreq, missing=missing)
	setMarkers(x, structure(c(x$markerdata, list(mm)), class="markerdata"))
}

setMarkers <- function(x, m, map=NULL, dat=NULL, freq=NULL, missing=0) {
	if(is.null(m)) 								markerdata_list = NULL
	else if(class(m) == "marker")				markerdata_list = structure(list(m), class="markerdata")
	else if(is.list(m) && all(unlist(lapply(m, class))=="marker")) markerdata_list = structure(m, class="markerdata")
	else if(class(m) == "markerdata")			markerdata_list = m
	else if ((n <- ncol(m <- as.matrix(m))) == 0)	markerdata_list = NULL
	else {
		is.merged = is.character(m[1,1]) && ("/" %in% strsplit(m[1,1],"")[[1]]); 
		if (!is.merged && n%%2 != 0) stop("Uneven number of marker allele columns")
		nMark <- ifelse(is.merged, n, n/2)
		
		if(is.null(map)) {
			markerdata_list = lapply(1:nMark, function(i) {
				if (is.merged)	mi = t(sapply(m[,i], function(markeri) strsplit(markeri,"/",fixed=T)[[1]]))
				else 			mi = m[, c(2*i-1, 2*i)]
				.createMarkerObject(mi, missing=missing, chrom=NA, pos=NA, name=NA)
		})}
		else {
			if(is.null(dat)) stop("Missing 'dat' file.") 
			mapinfo = .readMap(map, dat, freq)
			
			markerdata_list = lapply(seq_along(mapinfo), function(i) {
				info = mapinfo[i]; name=names(info); info=info[[1]]
				if(is.na(name)) return(NULL)
				
				if (is.merged)	mi = t(sapply(m[,i], function(markeri) strsplit(markeri,"/",fixed=T)[[1]]))
				else 			mi = m[, c(2*i-1, 2*i)]
				
				.createMarkerObject(mi, alleles=info$alleles, afreq=info$afreq, missing=missing, chrom=info$chrom, pos=info$pos, name=name) 
			})}
		markerdata_list[sapply(markerdata_list, is.null)] = NULL
		class(markerdata_list) = "markerdata"
	}
	x$nMark = length(markerdata_list) 
	x$markerdata = markerdata_list
	x$available = if (x$nMark>0) x$orig.ids[rowSums(do.call(cbind, markerdata_list)) > 0] else numeric(0)
	x
}

.readMap = function(map, dat, freq) {
	stopifnot(!is.null(map), !is.null(dat))
	if(is.character(map) && length(map)==1) {
		rawmap = read.table(map, as.is=TRUE, header=FALSE)
		if(is.character(rawmap[[1]]) && !any(1:9 %in% strsplit(rawmap[[1]][1],"")[[1]])) rawmap = rawmap[-1,] #If no number occurs in first entry, first row is assumed to be header.
	}
	else  rawmap = map 
	rawmap[[1]][rawmap[[1]]=="X"] = 23
	map1 = data.frame(CHR = as.numeric(rawmap[,1]), MARKER = as.character(rawmap[, 2]), POS = as.numeric(rawmap[, 3]),stringsAsFactors=FALSE) 

	if(is.character(dat) && length(map)==1)
		rawdat = read.table(dat, as.is=TRUE, header=FALSE)
	else  rawdat = dat
	datnames = as.character(rawdat[rawdat[, 1] == "M", 2])  #names of all markers in map
	
	freqlist=list()
	if(is.character(freq) && length(freq)==1) {
		rawfreq = read.table(freq, as.is=T, fill=T, header=FALSE) 
		if(rawfreq[2,1]=="F") {
			rawfreq = apply(rawfreq, 1, function(r) as.character(r[!is.na(r) & nzchar(r)]))
			for(i in 2*seq_len(length(rawfreq)/2))
				afreq = as.numeric(rawfreq[[i]][-1])
				freqlist[[rawfreq[[i-1]][2]]] = list(alleles = seq_along(afreq), afreq = afreq)
		}
		else { #if rawfreq[2,1]=="A" ... i.e. long format! See MERLIN tutorial about input files.
			Mrow = which(rawfreq[1]=="M")
			lastA = c(Mrow[-1] - 1, nrow(rawfreq))
			for(i in seq_along(Mrow)) {
				name = rawfreq[Mrow[i], 2]
				freqrows = (Mrow[i]+1):lastA[i]
				freqlist[[name]] = list(alleles = rawfreq[freqrows, 2], afreq = as.numeric(rawfreq[freqrows, 3]))
			}
		}
	}
	else if(is.list(freq))	freqlist = freq
	
	mapinfo = lapply(1:nrow(map1), function(i) list(chrom=map1$CHR[i], pos=map1$POS[i], alleles=freqlist[[map1[i,2]]]$alleles, afreq=freqlist[[map1[i,2]]]$afreq))
	names(mapinfo) = map1$MARKER
	
	Mmatch = match(datnames, map1$MARKER, nomatch=NA)
	if(any(NAs <- is.na(Mmatch)))
		cat("Deleting the following marker(s), which are not found in the map file:\n", paste(datnames[NAs], collapse="\n"), "\n", sep="")

	mapinfo[Mmatch]
}


.createMarkerObject = function(matr, missing, alleles=NULL, afreq=NULL, chrom=NA, pos=NA, name=NA) {
	if(is.null(alleles)) {
		vec = as.vector(matr)
		alleles = unique.default(vec[vec != missing])
		if (length(alleles)==0) alleles = 1
	}
	if(all(alleles %in% 0:99)) alleles = as.numeric(alleles)
	all_ord = order(alleles)
	alleles = as.character(alleles[all_ord])
	nalleles = length(alleles)
	if(is.null(afreq)) afreq = rep.int(1, nalleles)/nalleles 
	else {
		if(length(afreq) != nalleles) stop("Number of alleles don't match length of frequency vector")
		if(round(sum(afreq),2) != 1) warning(paste("Allele frequencies for marker", name," do not sum to 1:", paste(afreq, collapse=" ")))
		afreq = afreq[all_ord]
	}
	m_obj = match(matr, alleles, nomatch=0)
	attributes(m_obj) = list(dim=dim(matr), alleles=alleles, missing=missing, nalleles=nalleles, afreq = afreq, chrom=chrom, pos=pos, name=name, class="marker") 
	m_obj
}

marker = function(x, ..., allelematrix, alleles=NULL, afreq=NULL, missing=0, chrom=NA, name=NA, pos=NA) {
	arglist = list(...)
	n = length(arglist)
	if(n == 0) {
		if(missing(allelematrix)) m = matrix(missing, ncol=2, nrow=x$nInd)
		else {
			m = as.matrix(allelematrix)
			stopifnot(nrow(m)==x$nInd, ncol(m)==2)
		}
	}
	if(n > 0) {
		if(!missing(allelematrix)) stop("Syntax error. See ?marker.")
		if(n == 1 && length(arglist[[1]]) > 2) stop("Syntax error. See ?marker.")
		if(n > 1 && n%%2 != 0) stop("Wrong number of arguments.")
		
		fill = {if(n==1) arglist[[1]] else missing}
		m = matrix(fill, ncol=2, nrow=x$nInd, byrow=TRUE)  #create marker matrix with all individuals equal.
		
		for(i in (seq_len(n/2)*2 - 1)) {   # odd numbers up to n-1.
			ids = arglist[[i]]
			geno = arglist[[i+1]]
			for(j in match(ids, x$orig.ids)) m[j, ] = geno
		}
	}
	
	.createMarkerObject(m, missing=missing, alleles=alleles, afreq=afreq, chrom=chrom, name=name, pos=pos)
}

removeMarkers <- function(x, markers) {
	m = x$markerdata
	m[markers] = NULL
	setMarkers(x, m)
}

.prettyMarkers = function(m, alleles=NULL, sep="", missing=NULL, singleCol=FALSE, sex) {
	if(is.null(m)) return(m) 
	if(is.matrix(m)) m = list(m)
	if((n <- length(m))==0) return(m)
	if(is.null(alleles)) alleles = lapply(m, attr, 'alleles')
	else {
		if(!is.atomic(alleles)) stop("The parameter 'alleles' must be NULL, or an atomic vector.")
		if(length(alleles) < max(unlist(lapply(m, attr, 'nalleles')))) stop("The indicated 'alleles' vector has too few alleles for some markers.")
		alleles = rep(list(alleles), n)
	}
	if(is.null(missing)) missing = unlist(lapply(m, attr, 'missing')) 
	else {
		if(!is.atomic(missing) || length(missing) != 1) stop("The parameter 'mising' must be NULL, or a numeric/character of length 1.") 
		missing = rep(missing, n)
	}
	mNames = unlist(lapply(m, attr, 'name')); mNames[is.na(mNames)]=""
	pretty_m = do.call(c, lapply(seq_len(n), function(i)  c(missing[i], alleles[[i]])[m[[i]]+1]))
	dim(pretty_m) = c(length(pretty_m)/(2*n), 2*n)
	if (singleCol) {
		al1 = pretty_m[, 2*seq_len(n)-1, drop=FALSE]
		al2 = pretty_m[, 2*seq_len(n), drop=FALSE]
		m.matrix = matrix(paste(al1, al2, sep=sep), ncol=n )
		chrom_X = unlist(lapply(m, function(mm) identical(23L, as.integer(attr(mm, 'chrom')))))
		if(any(chrom_X)) {
			males = (sex==1)
			if (!all(hh <- al1[males, chrom_X] == al2[males, chrom_X])) warning("Male heterozygosity at X-linked marker detected.")
			m.matrix[males, chrom_X] = al1[males, chrom_X]
		}
		colnames(m.matrix) = mNames
		return(m.matrix)
	}
	else {
		nam = rep(mNames, each=2); nam[nzchar(nam)] = paste(nam[nzchar(nam)], 1:2, sep="_")
		colnames(pretty_m) = nam
		return(pretty_m)
	}
}

	
modifyMarker <- function(x, marker, ids, genotype, alleles, afreq, chrom, name, pos) {
	if(class(marker)=="marker") {
		if(nrow(marker) != x$nInd) stop("Wrong dimensions of marker matrix.")
		m = marker
	}
	else {
		if (!is.numeric(marker) || length(marker) != 1) stop("Argument 'marker' must be a single integer or an object of class 'marker'.")
		if (marker > x$nMark) stop("Indicated marker does not exist.")
		m = x$markerdata[[marker]]
	}
	mis = attr(m, 'missing')
	
	if(!missing(alleles)) {
		stopifnot(is.atomic(alleles), is.numeric(alleles) || is.character(alleles))
		#if(is.numeric(alleles)) alleles = as.character(alleles)
		if(attr(m, 'missing') %in% alleles) stop("The 'missing allele' character cannot be one of the alleles.")
		lena = length(alleles)
		if(lena == attr(m, 'nalleles'))
			attr(m, 'alleles') = alleles
		else {
			num_als = unique.default(as.vector(m[m!=0]))
			if(lena < length(num_als)) stop("Too few alleles.")
			pm = matrix(c(0, attr(m, 'alleles'))[m + 1], ncol=2)
			m = .createMarkerObject(pm, missing=0, alleles=alleles, afreq=rep(1, lena)/lena, chrom=attr(m, 'chrom'), pos=attr(m, 'pos'), name=attr(m, 'name'))
		}
	}
	if(!missing(afreq)) {
		if(round(sum(afreq),2)!=1) stop("The allele frequencies don't sum to 1.")
		if(length(afreq) != attr(m, 'nalleles')) stop("The length of allele frequency vector doesn't equal the number of alleles.")
		if(is.null(names(afreq))) 
			attr(m, 'afreq') = afreq
		else
			if (all(names(afreq) %in% attr(m, 'alleles'))) attr(m, 'afreq') = afreq[attr(m, 'alleles')] else stop("The names of the frequency vector don't match the allele names.")
	}
	
	changegeno = sum(!missing(ids), !missing(genotype))
	if(changegeno == 1) stop("The parameters 'ids' and 'genotype' must either both be NULL or both non-NULL.")
	if(changegeno == 2) {
		if(!is.atomic(genotype) || length(genotype)>2) stop("The 'genotype' parameter must be a numeric or character vector of length 1 or 2.")
		pm = .prettyMarkers(list(m))
		for (i in .internalID(x, ids)) pm[i, ] = genotype
		if(!all(as.vector(pm) %in% c(attr(m, 'alleles'), mis))) stop("Unknown allele(s). Please specify allele names using the 'alleles' argument.")
		attribs = attributes(m)
		m = match(pm, attr(m, 'alleles'), nomatch=0)
		attributes(m) = attribs
	}
	if(!missing(chrom)) attr(m, 'chrom') = chrom
	if(!missing(name)) attr(m, 'name') = name
	if(!missing(pos)) attr(m, 'pos') = pos
	if(class(marker)=="marker") return(m)
	else {
		x$markerdata[[marker]] = m
		x
	}
}

.getMarkers=function(x, markernames=NULL, chrom=NULL, startpos=NULL, endpos=NULL) {
	dat = data.frame(CHR = sapply(x$markerdata, attr, 'chrom'), POS = sapply(x$markerdata, attr, 'pos'), NAME = sapply(x$markerdata, attr, 'name'), stringsAsFactors=F)
	test = !logical(nrow(dat)) 
	if(!is.null(markernames)) test = test & (dat$NAME %in% markernames)
	if(!is.null(chrom)) test = test & (dat$CHR %in% chrom)
	if(!is.null(startpos)) test = test & (dat$POS > startpos)
	if(!is.null(endpos)) test = test & (dat$POS < endpos)
	which(test)
}

.swapMarkerData=function(x, ids){ #Laget for Thore, eksporteres ikke.
	stopifnot(length(ids)==2)
	ids = .internalID(x, ids)
	y = .as.annotated.matrix(x)
	y[ids, -(1:6)] = y[ids[2:1], -(1:6)]
	.restore.linkdat(y)
}

.modifyMarkerData=function(x, ids, new.alleles){
	ids = .internalID(x, ids)
	y = .as.annotated.matrix(x)
	y[ids, -(1:6)] = new.alleles
	.restore.linkdat(y)
}


.diallel2geno <- function(marker) { #marker a numerical nInd * 2 matrix   
	#Coding genotypes as single integer: 00 -> 0, 11 -> 1, 22 -> 2, 12/21 -> 3, 01/10 -> 4, 02/20 -> 5.
	#Each pair of alleles is seen as an integer written in base 3, and this integer is permuted to fit with the above code.
	c(0,4,5,4,1,3,5,3,2)[colSums(c(3,1) * t(marker)) + 1]
}

.geno2diallel <- function(codedgenos) { #input: matrix of single-numerical genotypes. 
	#Ouput: Matrix with twice the number of columns, decoded as 1 -> 1 1, 2 -> 1 2, 3 -> 2 2, 4 -> 1 0, 5 -> 2 0
	decode = sapply(t(codedgenos+1),function(i) switch(i, c(0,0), c(1,1), c(2,2), c(1,2), c(1,0), c(2,0)))
	dim(decode) = c(2*ncol(codedgenos), nrow(codedgenos))
	t(decode)
}


mendelianCheck <- function(x, verbose=TRUE) {

	trioCheckFast = function(fa, mo, of) {
		fa_odd = fa[odd]; fa_even = fa[even]; mo_odd = mo[odd]; mo_even = mo[even]
		of_odd = of[odd]; of_even = of[even];
		fa0 = (fa_odd == 0 | fa_even == 0)
		mo0 = (mo_odd == 0 | mo_even == 0)
		of_odd0 = (of_odd == 0); of_even0 = (of_even == 0)
		ff1 = (fa0 | of_odd0 | of_odd == fa_odd | of_odd == fa_even)
		ff2 = (fa0 | of_even0 | of_even == fa_odd | of_even == fa_even)
		mm1 = (mo0 | of_odd0 | of_odd == mo_odd | of_odd == mo_even)
		mm2 = (mo0 | of_even0 | of_even == mo_odd | of_even == mo_even)
		(ff1 & mm2) | (ff2 & mm1)
	}
	
	maleXHomoz = function(of) 	of[odd] == of[even]
	
	maleXCheck = function(mo, of)	{
		mo_odd = mo[odd]; mo_even = mo[even]
		of = of[odd]
		mo0 = (mo_odd == 0 | mo_even == 0)
		of == 0 | mo0 | of == mo_odd | of == mo_even
	}
	
	nuclearAllelCheck = function(parents, offs) { #both matrices with 2*N columns
		unlist(lapply(even, function(i) {
			offs_als = unique.default(offs[,(i-1):i])
			offs_count = length(offs_als[offs_als > 0])
	
			par_als = parents[,(i-1):i]
			iszero = (par_als == 0); nzero = sum(iszero)
			if(nzero==0) 	parents_count = length(unique.default(par_als))
			else { 
				par_als = unique.default(par_als[!iszero])
				parents_count = length(.myintersect(par_als, offs_als)) + nzero
			}
			offs_count <= parents_count
		}))
	}	
	
	ped = x$pedigree
	parents = unique(ped[, 2:3]) 
	parents = parents[ -match(0, parents[,1]), , drop=FALSE]
	subnucs = lapply(nrow(parents):1, function(i) { 
		par = parents[i,];  
		c(fa=par[[1]], mo=par[[2]], offs=as.vector(ped[,1])[ which(ped[,2]==par[[1]] & ped[,3]==par[[2]], useNames=FALSE) ]) 
	})
	
	chroms = ifelse(all(unlist(lapply(x$markerdata, function(mm) identical(23L, as.integer(attr(mm, 'chrom')))))), 'X', 'AUTOSOMAL')
	mdat = do.call(cbind, x$markerdata)
	Nseq = seq_len(length(x$markerdata)); even = 2* Nseq; odd = even - 1 
	
	error = numeric(0)
	switch(chroms, 
	AUTOSOMAL = {
		for (sub in subnucs) {
			fa = mdat[sub[[1]], ]; mo = mdat[sub[[2]], ]; offs = sub[-(1:2)]
			for (of in offs) {
				new_errors = Nseq[!trioCheckFast(fa, mo, mdat[of, ])]
				if(verbose) for (i in new_errors) cat("Marker ", i, ", individual ", x$orig.ids[of], ": Alleles not compatible with parents.\n", sep="")
				error = c(error, new_errors)
			}
			if (length(offs) > 1) {
				new_errors = Nseq[!nuclearAllelCheck(parents = mdat[sub[1:2], ], offs = mdat[sub[-(1:2)], ])]
				if(verbose) for (i in new_errors)  cat("Marker ", i, ": Offspring of ", x$orig.ids[sub[1]], " and ", x$orig.ids[sub[2]], " have too many different alleles.\n", sep="")
				error = c(error, new_errors)
			}
		}
	}, 
	X = {
	sex = x$pedigree[, 'SEX']
	for (sub in subnucs) {
			fa = mdat[sub[[1]], ]; mo = mdat[sub[[2]], ]; offs = sub[-(1:2)]
			for (of in offs) {
				ofdat = mdat[of, ]
				if(sex[of]==1) {
					new_errors = Nseq[!maleXHomoz(ofdat)]
					if(verbose) for (i in new_errors) cat("Marker ", i, ", individual ", x$orig.ids[of], ": Male heterozygosity not compatible with X-linked model.\n", sep="")
					error = c(error, new_errors)
					
					new_errors = Nseq[!maleXCheck(mo, ofdat)]
					if(verbose) for (i in new_errors) cat("Marker ", i, ", individual ", x$orig.ids[of], ": Allele not compatible with mother.\n", sep="")
					error = c(error, new_errors)
				}
				else {
					new_errors = Nseq[!trioCheckFast(fa, mo, ofdat)]
					if(verbose) for (i in new_errors) cat("Marker ", i, ", individual ", x$orig.ids[of], ": Alleles not compatible with parents.\n", sep="")
					error = c(error, new_errors)
				}
			}
			if (length(offs) > 2) {
				new_errors = Nseq[!nuclearAllelCheck(parents = mdat[sub[1:2], ], offs = mdat[offs, ])] #TODO: fix this for X-linked
				if(verbose) for (i in new_errors)  cat("Marker ", i, ": Offspring of ", x$orig.ids[sub[1]], " and ", x$orig.ids[sub[2]], " have too many different alleles.\n", sep="")
				error = c(error, new_errors)
			}
		}
	})
	sort.int(unique.default(error))
}
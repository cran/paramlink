linkdat <-
function(ped, model=NULL, map=NULL, dat=NULL, freq=NULL, verbose=TRUE, missing=0, ...) {  

	subnucs <- function(ped) {  	# output: peeling order of nuclear subfamilies. Format for each nuc: c(pivot,father,mother,offsp1,..), where pivot=0 for the last nuc.
		parents = unique(ped[, 2:3])
		parents = parents[ -match(0, parents[,1]), , drop=FALSE]
		list1 = lapply(nrow(parents):1, function(i) { 
			par = parents[i,]; 
			list(father = par[[1]], mother = par[[2]], offspring = as.vector(ped[,1])[ which(ped[,2]==par[[1]] & ped[,3]==par[[2]], useNames=FALSE) ]) } )  #listing all nucs
		res = list(); i = 1; k = 1
		while(length(list1) > 1)	{
			if (i > length(list1)) 	return(FALSE)
			sub = list1[[i]]
			subvec = unlist(sub)
			links = subvec[subvec %in% unlist(list1[-i])]
			if (length(links)==1) 	{
				res[[k]] <- c(sub, list(pivot = as.numeric(links), pivtype = match(links, c(sub[['father']], sub[['mother']]), nomatch=3)))
				list1 <- list1[-i]
				k <- k+1;	i <- 1
			} else i <- i+1
		}	
		res[[k]] <- c(list1[[1]], list(pivot = 0, pivtype = 0)) #final nuclear
		res
	}
	
	if (class(ped)=="linkdat") 
		return(ped)	
	if(is.character(ped) && length(ped)==1)
		ped = read.table(ped, as.is=TRUE, ...)
	if(nrow(ped) < 3) stop("The pedigree must have at least 3 individuals.")
	
	if(length(famids <- unique(ped[, 1])) < nrow(ped)) {
		if(length(famids) > 1) 
			return(lapply(famids, function(fam)  linkdat(ped[ ped[,1]==fam, ], model=model, map=map, dat=dat, freq=freq, verbose=verbose, missing=missing)))
		else {
			famid = as.vector(ped[1,1])
			ped=ped[,-1]
		}
	}
	else famid=1
	if (verbose) cat("Family name: ", famid, ".\n", sep="")
	if(ncol(ped) < 5) stop("Too few columns: ID, FID, MID, SEX and AFF are mandatory.")
	
	pedigree = as.matrix(ped[, 1:5], rownames.force=FALSE)
	colnames(pedigree) = c('ID', 'FID', 'MID', 'SEX', 'AFF')
	.checkped(pedigree)
	
	orig.ids = ped[,1]
	nInd = nrow(pedigree)
	pedigree = relabel(pedigree, new=1:nInd)
	if (verbose) cat("Pedigree read,", nInd, "individuals.\n")
	
	founders = as.integer(which(pedigree[,'FID'] == 0))
	nonfounders = as.integer(which(pedigree[,'FID'] > 0))
	
	#---peeling order of nuclear subfamilies---
	subnucs = subnucs(pedigree)
	hasLoops = is.logical(subnucs) && !subnucs
	if(verbose) 
		if(hasLoops) cat("Loop(s) detected.\n")  
		else if(is.list(subnucs)) cat(ant<-length(subnucs), "maximal nuclear", ifelse(ant==1,"subfamily.\n", "subfamilies. Peeling order set.\n"))
	
	#---creation of linkdat object---
	obj = structure(list(pedigree=pedigree, famid=famid, orig.ids=orig.ids, nInd=nInd, founders=founders, 
			nonfounders=nonfounders, hasLoops=hasLoops, subnucs=subnucs), class="linkdat")

	#---adding markers---
	markercols = seq_len(ncol(ped)-5)+5
	obj = setMarkers(obj, as.matrix(ped[, markercols], rownames.force=FALSE), map=map, dat=dat, freq=freq, missing=missing)
	if (verbose) if(obj$nMark==1) cat("1 marker read.\n") else cat(obj$nMark, "markers read.\n")

	#---adding availability vector---
		#this is now done inside setMarkers()

	#----adding model----
 	if (!is.null(model)) obj = setModel(obj, model=model)

	#----adding map----
 	#if (!is.null(map)) obj = setMap(obj, map, dat, verbose=verbose)

	invisible(obj)
}

.checkped <- function(p) { #p a numeric matrix with 5 columns
	if (!is.numeric(p)) stop("Pedigree columns must be numeric.")
	ID = p[,'ID']; FID = p[,'FID']; MID = p[,'MID']; SEX = p[,'SEX']; AFF=p[,'AFF']
	if (all(c(FID, MID)==0)) stop("Pedigree is not connected.")

	quick.check <- (!any(duplicated(ID)) && all(ID!=FID) && all(ID!=MID) && all(c(FID, MID) %in% c(0,ID)) && 
					all(SEX[match(FID[FID!=0], ID)] == 1) && all(SEX[match(MID[MID!=0], ID)] == 2) && 
					all(SEX %in% 1:2) && all(AFF %in% 0:2) && all((FID>0)==(MID>0)))
	
	if (quick.check)  #if all tests are passed
		return()
	else {
	  for (i in seq_along(ID)) {
		if (i>1 && ID[i] %in% ID[1:(i-1)]) cat("Individual ", ID[i],": ID not unique.\n")
		if (!FID[i] %in% c(0,ID)) cat("Individual ", ID[i],": Father's ID does not appear in ID column.\n")
		if (!MID[i] %in% c(0,ID)) cat("Individual ", ID[i],": Mother's ID does not appear in ID column.\n")
		if (!SEX[i] %in% 1:2) cat("Individual ", ID[i],": SEX must be either 1 (male) or 2 (female).\n")
		if (!AFF[i] %in% 0:2) cat("Individual ", ID[i],": Affection status must be either 0 (unknown), 1 (non-affected) or 2 (affected).\n")
		if (FID[i] != 0 && SEX[match(FID[i], ID)] != 1) cat("Individual ", ID[i],": Father is not male.\n")
		if (MID[i] != 0 && SEX[match(MID[i], ID)] != 2) cat("Individual ", ID[i],": Mother is not female.\n")
		if ((FID[i]>0) != (MID[i]>0)) cat("Individual ", ID[i],": Only one parent in the pedigree is not allowed. Either both parents or none must be specified.\n")
		if (ID[i]==FID[i]) cat("Individual ", ID[i]," is his own father. Paramlink does not accept this sort of behaviour.\n")
		if (ID[i]==MID[i]) cat("Individual ", ID[i]," is her own mother. Paramlink does not accept this sort of behaviour.\n")
	  }
	  stop("Pedigree errors detected.")
	}
}

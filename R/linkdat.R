linkdat <-
function(x, model=NULL, verbose=TRUE, missing=0) {  

	subnucs <- function(ped) {  	# output: peeling order of nuclear subfamilies. Format for each nuc: c(pivot,father,mother,offsp1,..), where pivot=0 for the last nuc.
		parents <- unique(ped[, 2:3]); parents = parents[ -match(0, parents[,1]), , drop=FALSE]
		list1 <- lapply(nrow(parents):1, function(i) { par=parents[i,]; c(fa = par[[1]], mo = par[[2]], offs = as.vector(ped[,1])[ which(ped[,2]==par[[1]] & ped[,3]==par[[2]], useNames=FALSE) ]) } )  #listing all nucs
		res=list(); i=1; k = 1
		while(length(list1) > 1)	{
			if (i > length(list1)) {cat("Loop(s) detected!\n"); return(NULL)} 
			link = ( (sub <- list1[[i]]) %in% unlist(list1[-i]))
			if (sum(link)==1) 	{
				res[[k]] <- c(pivot=sub[[which(link)]], sub)
				list1 <- list1[-i]
				k <- k+1;	i <- 1
			} else i <- i+1
		}	
		res[[k]] <- c(pivot=0, list1[[1]]) #final nuclear
		res
	}
  
	if (!class(x) %in% c("linkdat", "matrix", "data.frame")) 
		stop("Input must be either a matrix, a data.frame or a 'linkdat' object.")
	if (class(x)=="linkdat") 
		return(x) #x = as.data.frame(x, famid=T) 
	
	if(nrow(x) < 3) stop("The pedigree must have at least 3 individuals.")
	if(x[1,1]==x[2,1]) {famid = x[1,1]; x=x[,-1]; if (verbose) cat("Interpreting first column as family identifier\n")} 
	else famid=1
	if(ncol(x) < 5) stop("Too few columns: ID, FID, MID, SEX and AFF are mandatory.")
	
	pedigree = as.matrix(x[, 1:5], rownames.force=FALSE)
	colnames(pedigree) = c('ID', 'FID', 'MID', 'SEX', 'AFF')
	.checkped(pedigree)
	
	orig.ids = x[,1]
	nInd = nrow(pedigree)
	pedigree = relabel(pedigree, new.labels=1:nInd)
	if (verbose) cat("Pedigree read,", nInd, "individuals\n")
	
	founders = as.integer(which(pedigree[,'FID'] == 0));	nonfounders = as.integer(which(pedigree[,'FID'] > 0))
	
	#---peeling order of nuclear subfamilies---
	subnucs = subnucs(pedigree)
	if (verbose && !is.null(subnucs)) cat(ant<-length(subnucs), "maximal nuclear", ifelse(ant==1,"subfamily.\n", "subfamilies. Peeling order set.\n"))
	
	#---creation of linkdat object---
	obj = structure(list(pedigree=pedigree, famid=famid, orig.ids=orig.ids, nInd=nInd, founders=founders, 
			nonfounders=nonfounders, subnucs=subnucs), class="linkdat")

	#---adding markers---
	markercols = seq_len(ncol(x)-5)+5
   	obj = setMarkers(obj, as.matrix(x[, markercols], rownames.force=FALSE), missing=missing)
    if (verbose) cat(obj$nMark, "markers read.\n")
	
	#---if marker data exists, set simulation vector (those that are genotyped with at least one marker, are set as TRUE)
	if (obj$nMark > 0) obj = setSim(obj, rowSums(obj$markerdata) > 0)
	
	#----adding model----
 	if (!is.null(model)) obj = setModel(obj, model=model)
	
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

.linkdat_new <-
function(file, pedigree, model=NULL, markerdata=NULL, map=NULL, sim=NULL, verbose=TRUE, missing=0) {  

	subnucs <- function(ped) {  	# output: peeling order of nuclear subfamilies. Format for each nuc: c(pivot,father,mother,offsp1,..), where pivot=0 for the last nuc.
		parents <- unique(ped[, 2:3]); parents = parents[ -match(0, parents[,1]), , drop=FALSE]
		list1 <- lapply(nrow(parents):1, function(i) { par=parents[i,]; c(fa = par[[1]], mo = par[[2]], offs = as.vector(ped[,1])[ which(ped[,2]==par[[1]] & ped[,3]==par[[2]], useNames=FALSE) ]) } )  #listing all nucs
		res=list(); i=1; k = 1
		while(length(list1) > 1)	{
			if (i > length(list1)) {warning("Loop detected, likelihood calculations will not work."); return(NULL)} 
			link = ( (sub <- list1[[i]]) %in% unlist(list1[-i]))
			if (sum(link)==1) 	{
				res[[k]] <- c(pivot=sub[[which(link)]], sub)
				list1 <- list1[-i]
				k <- k+1;	i <- 1
			} else i <- i+1
		}	
		res[[k]] <- c(pivot=0, list1[[1]]) #final nuclear
		res
	}
  
	if (!class(x) %in% c("linkdat", "matrix", "data.frame")) 
		stop("Input must be either a matrix, a data.frame or a 'linkdat' object.")
	if (class(x)=="linkdat") 
		return(x) #x = as.data.frame(x, famid=T) 
	
	if(x[1,1]==x[2,1]) {famid = x[1,1]; x=x[,-1]; if (verbose) cat("Interpreting first column as family identifier\n")} 
	else famid=1

	pedigree = as.matrix(x[, 1:5], rownames.force=FALSE)
	colnames(pedigree) = c('ID', 'FID', 'MID', 'SEX', 'AFF')
	.checkped(pedigree)
	
	orig.ids = x[,1]
	nInd <- nrow(pedigree)
	pedigree <- relabel(pedigree, new.labels=1:nInd)
	if (verbose) cat("Pedigree read,", nInd, "individuals\n")
	
	founders <- as.integer(which(pedigree[,'FID'] == 0));	nonfounders <- as.integer(which(pedigree[,'FID'] > 0))

	# #---simulation column---
	# simcol <- match('SIM', colnames(x))
	# if (!is.na(simcol)) {
		# sim <- as.integer(x[, simcol])
		# if (verbose) cat("Simulation column read\n")
	# } else {
		# sim=NULL 
		# if (verbose) cat("No simulation column found\n")
	# }
	
	#---peeling order of nuclear subfamilies---
	subnucs = subnucs(pedigree)
	if (verbose) cat("Peeling order set,", ant<-length(subnucs), "maximal nuclear", ifelse(ant==1,"subfamily.","subfamilies."),"\n")
	
	#---creation of linkdat object--
	obj <- structure(list(pedigree=pedigree, famid=famid, orig.ids=orig.ids, nInd=nInd, founders=founders, 
			nonfounders=nonfounders, subnucs=subnucs), class="linkdat")

	#---adding markers
	markercols <- seq_len(ncol(x)-5)+5
   	obj = setMarkers(obj, as.matrix(x[, markercols], rownames.force=FALSE), missing=missing)
    if (verbose) cat(obj$nMark, "markers read.\n")
	
	#----adding model
 	if (!is.null(model)) obj <- setModel(obj, model=model)
	
	invisible(obj)
}

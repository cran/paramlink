linkdat <-
function(x, model=NULL, verbose=TRUE, missing=0) {  

	subnucs <- function(ped) {  	# output: peeling order of nuclear subfamilies. Format for each nuc: c(pivot,father,mother,offsp1,..), where pivot=0 for the last nuc.
		parents <- unique(ped[, 2:3]); parents = parents[ -match(0, parents[,1]), , drop=FALSE]
		list1 <- lapply(nrow(parents):1, function(i) { par=parents[i,]; c(fa = par[[1]], mo = par[[2]], offs = as.vector(ped[,1])[ which(ped[,2]==par[[1]] & ped[,3]==par[[2]], useNames=FALSE) ]) } )
		res=list(); i=1; k = 1
		while(length(list1) > 1)	{
			if (i > length(list1)) {warning("Loop detected, likelihood calculations will not work."); return(NULL)} 
			link = ( (sub<-list1[[i]]) %in% unlist(list1[-i]))
			if (sum(link)==1) 	{
				res[[k]] <- c(pivot=sub[[which(link)]], sub)
				list1 <- list1[-i]
				k <- k+1;	i <- 1
			} else i <- i+1
		}	
		res[[k]] <- c(pivot=0, list1[[1]]) #final nuclear
		res
	}
  
	renumerate <- function(df, new_labels) {
		orig.ids <- df[, 'ID']

		#identify and remove parents outside of the pedigree (usually none, except when using subset.linkdat().)
		parents.out <- apply(df[,c('FID','MID')], 1, function(par) any(!par %in% c(0, orig.ids)))
		df[ parents.out, c('FID','MID')] <- c(0,0)

		df[,'ID'] <- new_labels
		
		#relabeling parents
		for (i in 1:nrow(df)) { 
			par=df[i, c('FID', 'MID')]
			if ( any( par > 0 ) ) df[i, c('FID', 'MID')] <- new_labels[match(par, orig.ids)]
		}
		if (verbose) cat("Relabeled individuals:\n", paste( paste(orig.ids, new_labels, sep=" -> "), "\n", sep=""))
		return(df)
	}
  
	if (!class(x) %in% c("linkdat", "matrix", "data.frame")) 
		stop("Input must be either a matrix, a data.frame or a 'linkdat' object.")
	if (class(x)=="linkdat") 
		x = as.data.frame(x) 
	
	pedigree = as.matrix(x[, 1:5], rownames.force=FALSE)
	.checkped(pedigree)
	nInd = nrow(pedigree)
	if (verbose) cat("Pedigree read,", nInd, "individuals\n")
	
	oldid <- x[, match('OLDID', colnames(x), nomatch=1)]
	
	if (!all( pedigree[,'ID'] == 1:nInd)) 
		pedigree <- renumerate(pedigree, new_labels=1:nInd)
	
	founders <- as.integer(which(pedigree[,'FID'] == 0));	nonfounders <- as.integer(which(pedigree[,'FID'] > 0))

	#---simulation column---
	simcol <- match('SIM', colnames(x))
	if (!is.na(simcol)) {
		sim <- as.integer(x[, simcol])
		if (verbose) cat("Simulation column read\n")
	} else {
		sim=NULL 
		if (verbose) cat("No simulation column found\n")
	}
	
	#---peeling order of nuclear subfamilies---
	subnucs = subnucs(pedigree)
	if (verbose) cat("Peeling order set,", ant<-length(subnucs), "maximal nuclear", ifelse(ant==1,"subfamily.","subfamilies."),"\n")
	
	#---creation of linkdat object--
	obj <- structure(list(pedigree=pedigree, oldid=oldid, nInd=nInd, founders=founders, 
			nonfounders=nonfounders, subnucs=subnucs, sim=sim), class="linkdat")

	#---adding markers
	markercols <- which(!colnames(x) %in% c('ID','FID','MID','SEX','AFF','OLDID','SIM'))
	obj = setMarkers(obj, as.matrix(x[, markercols], rownames.force=FALSE))
	if (verbose) cat(obj$nMark, "markers read.\n")
		
	#----adding model
 	if (!is.null(model)) obj <- setModel(obj, model=model)
	
	invisible(obj)
}

.checkped <- function(p) { #p a numeric matrix with 5 columns
	#p = as.matrix(x)
	if (!is.numeric(p)) stop("Pedigree columns must be numeric.")
	if (any(colnames(p) != c('ID','FID','MID','SEX','AFF'))) stop("Pedigree columns must be named 'ID', 'FID', 'MID', 'SEX', 'AFF'.")
	ID = p[,'ID']; FID = p[,'FID']; MID = p[,'MID']; SEX = p[,'SEX']; AFF=p[,'AFF']
	quick.check <- (!any(duplicated(ID)) && all(ID!=FID) && all(ID!=MID) && all(c(FID, MID) %in% c(0,ID)) && 
					all(SEX[match(FID[FID!=0], ID)] == 1) && all(SEX[match(MID[MID!=0], ID)] == 2) && 
					all(SEX %in% 1:2) && all(AFF %in% 0:2))
	
	if (quick.check)  #if all tests are passed
		return()
	else for (i in seq_along(ID)) {
		if (i>1 && i %in% ID[1:(i-1)]) stop("\nIndividual ", ID[i],": ID not unique.")
		if (!FID[i] %in% c(0,ID)) stop("\nIndividual ", ID[i],": Father's ID does not appear in ID column.")
		if (!MID[i] %in% c(0,ID)) stop("\nIndividual ", ID[i],": Mother's ID does not appear in ID column.")
		if (!SEX[i] %in% 1:2) stop("\nIndividual ", ID[i],": SEX must be either 1 or 2.")
		if (!AFF[i] %in% 0:2) stop("\nIndividual ", ID[i],": Affection status must be either 0, 1 or 2.")
		if (FID[i] != 0 && SEX[match(FID[i], ID)] != 1) stop("\nIndividual ", ID[i],": Father is not male.")
		if (MID[i] != 0 && SEX[match(MID[i], ID)] != 2) stop("\nIndividual ", ID[i],": Mother is not female.")
		if (ID[i]==FID[i]) stop("\nIndividual ", ID[i],": Illigal ID of father.")
		if (ID[i]==MID[i]) stop("\nIndividual ", ID[i],": Illigal ID of mother.")
	}	
}

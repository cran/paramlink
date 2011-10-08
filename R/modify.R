relabel <- function(x, new.labels) {
	if(islinkdat <- (class(x)=="linkdat")) ped = as.data.frame(x, famid=T) else ped = x
	stopifnot(is.numeric(new.labels), length(new.labels)==nrow(ped), !0 %in% new.labels)
	orig.ids = ped[, 'ID']
	ped[, 'ID'] = new.labels
	nonzero_par = ped[, c('FID','MID')] > 0
	ped[, c('FID','MID')][nonzero_par] <- new.labels[match(ped[, c('FID','MID')][nonzero_par], orig.ids)] #relabeling parents
	if(islinkdat) {
		mis = ifelse(is.null(x$markerdata), 0, attr(x$markerdata, "missing"))
		ped = linkdat(ped, model=x$model, verbose=FALSE, missing=mis)
		ped = setSim(ped, x$sim)
	}
	ped
}

setSim = function(x, simstatus) {
	if (is.null(simstatus)) {x$sim = NULL; return(x)}
	if(length(simstatus)==1) simstatus = rep(simstatus, x$nInd) else stopifnot(length(simstatus) == x$nInd)    
	if(is.logical(simstatus)) {sim = rep.int(0, x$nInd); sim[simstatus] = 2}	else sim = simstatus
	if(!all(sim %in% c(0,2)))  stop("The argument 'simstatus' must be either a logical vector, or a numeric consisting of 0's and 2's.")
	x$sim = sim
	x
}

swapSim = function(x, ids) {
	if(is.null(x$sim)) stop("Simulation vector not set. Use 'setSim()' first.")
	ids = .internalID(x, ids)
	x$sim[ids] = 2 - x$sim[ids]
	x
}
	
swapSex <- function(x, ids) {
	ids = .internalID(x, ids)
	ids.spouses = unique(unlist(lapply(ids, spouses, x=x, original.id = FALSE)))
	if (!all(ids.spouses %in% ids)) {
		cat("Changing sex of the following spouses as well:", paste(x$orig.ids[setdiff(ids.spouses, ids)], collapse=", "), "\n")
		return(swapSex(x, x$orig.ids[union(ids, ids.spouses)]))
	}	
	df=as.data.frame(x, famid=T)
	df[ids, 'SEX'] <- (3 - df[ids, 'SEX'])
	offs = x$pedigree[,"FID"] %in% ids
	df[offs, c('FID', 'MID')] <- df[offs, c('MID', 'FID')]
	
	mis = ifelse(is.null(x$markerdata), 0, attr(x$markerdata, "missing"))
	newx = linkdat(df, model=x$model, verbose=FALSE, missing=mis)
	setSim(newx, x$sim)
}

swapAff <- function(x, ids, newval=NULL) {
	df = as.data.frame(x, famid=T)
	ids = .internalID(x, ids)
	if (is.null(newval))  newval <- (3 - df[ids, 'AFF'])
	
	df[ids, 'AFF'] <- newval
	
	mis = ifelse(is.null(x$markerdata), 0, attr(x$markerdata, "missing"))
	newx = linkdat(df, model=x$model, verbose=FALSE, missing=mis)
	setSim(newx, x$sim)
}


addOffspring <- function(x, father, mother, noffs, ids, sex=1, aff=1) {
	p = relabel(x$pedigree, x$orig.ids)
	taken <- oldids <- p[,'ID']
	if(!missing(father)) taken = c(taken, father)
	if(!missing(mother)) taken = c(taken, mother)
	if(!missing(ids)) taken = c(taken, ids)
	max_id = max(taken)
	
	if(missing(father) && missing(mother)) stop("At least one parent must be an existing pedigree member.")
	if(missing(father)) 	father <- max_id <- max_id + 1 
	if(missing(mother)) 	mother <- max_id <- max_id + 1 
	if(any(!is.numeric(father), length(father)!=1)) 	stop("Argument 'father' must be a single integer.")
	if(any(!is.numeric(mother), length(mother)!=1)) 	stop("Argument 'mother' must be a single integer.")
	if(!any(c(father,mother) %in% oldids))	stop("At least one parent must be an existing pedigree member.")

	if (missing(noffs) && missing(ids)) stop("Number of offspring not indicated.")
	if (missing(noffs)) noffs = length(ids)
	if (missing(ids)) ids = (max_id+1):(max_id+noffs)
	if (length(ids)!=noffs) stop("Length of 'id' vector must equal number of offspring.")
	if (any(ids %in% oldids))	stop(paste("Individual(s)", ids[ids %in% oldids], "already exist(s)."))

	if(!father %in% oldids) 	p = rbind(p, c(father, 0, 0, 1, 1))
	if(!mother %in% oldids) 	p = rbind(p, c(mother, 0, 0, 2, 1))
	
	p = rbind(p, cbind(ids, father, mother, sex, aff))
	new_x = linkdat(p, model=x$model, verbose=FALSE)
	
	if (!is.null(m <- x$markerdata)) {
		mis = attr(m, "missing")
		new_m = rbind(.prettyMarkers(m), matrix(mis, nrow=new_x$nInd-x$nInd, ncol=ncol(m)))
		new_x = setMarkers(new_x, new_m, missing=mis)
	}
	new_x
}

addParents <- function(x, id, father, mother) {
	if(length(id)>1) stop("Only one individual at the time, please")
	if(id %in% x$orig.ids[x$nonfounders]) stop(paste("Individual", id, "already has parents in the pedigree")) 
	
	p = relabel(x$pedigree, x$orig.ids)
	n = max(p[,'ID'])
	if (missing(father)) father <- n <- n+1	
	if (missing(mother)) mother = n + 1

	p[id.row_in <- match(id, p[,'ID']), 2:3] <- c(father, mother)
	
	if(!father %in% p[,'ID']) {
		cat("Father: Creating new individual with ID", father, "\n")
		p = rbind(p, c(father,0,0,1,1))[append(1:nrow(p), nrow(p)+1, after=id.row_in-1), ] #insert row
	}
	if(!mother %in% p[,'ID']) {
		cat("Mother: Creating new individual with ID", mother, "\n")
		p = rbind(p, c(mother,0,0,2,1))[append(1:nrow(p), nrow(p)+1, after=match(id, p[,'ID'])-1), ] #insert row
	}

	new_x = linkdat(p, model=x$model, verbose=FALSE)
	
	if (!is.null(m <- x$markerdata)) {
		new_m = .prettyMarkers(m); mis = attr(m, "missing"); nro_m = nrow(new_m)
		if (0 < (added <- nrow(p) - nro_m)) { #if any individuals has been added
			new_m = rbind(new_m, matrix(mis, nrow=added, ncol=ncol(m)))[append(1:nro_m, nro_m + 1:added, id.row_in-1), ] #insert row(s)
		}
		new_x = setMarkers(new_x, new_m, missing=mis)
	}
	new_x
}

removeIndiv <- function(x, ids) { #removes (one by one) individuals 'ids' and all their descendants. Spouse-founders are removed as well.
	if(length(ids)==0) return(x)
	if(any(!ids %in% x$orig.ids)) stop(paste("Some of the indicated individuals do not exist:", paste(ids[!ids %in% x$orig.ids], collapse=", ")))  
	id = ids[[1]]
	internal_id = .internalID(x, id)

	#founders without children after 'id' and 'desc' indivs are removed. The redundancy here does not matter.
	desc = descendants(x, internal_id, original.id=F)
	leftover.spouses = setdiff(x$founders, c(internal_id, as.numeric(x$pedigree[ -c(internal_id, desc) , c('FID','MID')])))   #last part: remaining parents
	remov = c(internal_id, desc, leftover.spouses)
	cat("Removing individual", id, "and all his/her descendants (including leftover spouses):", paste(x$orig.ids[c(desc, leftover.spouses)], collapse=", "), "\n")
	new_df = as.data.frame(x, famid=T)[-remov, ] 

	mis = ifelse(is.null(x$markerdata), 0, attr(x$markerdata, "missing"))
	new_x = linkdat(new_df, model=x$model, missing=mis, verbose=FALSE)
	new_x = setSim(new_x, x$sim[-remov])
	removeIndiv(new_x, ids=ids[!ids %in% x$orig.ids[remov]])
}
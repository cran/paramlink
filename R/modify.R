modifyPedigree <- function(x, id, attrib, newval=NULL) {
	attrib=toupper(attrib)
	if (!attrib %in% c('AFF','SEX')) 
		stop("The 'attrib' argument must be either 'SEX' or 'AFF'.")
	if (attrib=='SEX' && length(spouses(x,id))>0) 
	stop("Sorry, gender change of individuals with offspring is not implemented at the moment.")
	df=as.data.frame(x)
	if (is.null(newval))
		newval <- (3 - df[id, attrib])
	else stopifnot(length(newval)==length(id), newval %in% c(0,1,2))

	df[id, attrib] <- newval
	mis = ifelse(is.null(x$markerdata), 0, attr(x$markerdata, "missing"))
	linkdat(df, model=x$model, verbose=FALSE, missing=mis)
}

addChildren <- function(x, father, mother, children=1, sex=1, aff=1) {
	if (!is.numeric(c(father, mother)) || length(c(father,mother))!=2) stop("Both arguments 'father' and 'mother' must be non-negative integers.")
	if (father==0 && mother==0) stop("At least one parent must be an existing pedigree member.")
	p = x$pedigree
	n = x$nInd
	if (father==0) 
		p = rbind(p, c(father <- (n+1), 0, 0, 1, 1))
	if (mother==0)
		p = rbind(p, c(mother <- (n+1), 0, 0, 2, 1))
	
	p = rbind(p, cbind(seq.int(nrow(p)+1, length=children), father, mother, sex, aff))
	new_x = linkdat(p, model=x$model, verbose=FALSE)
	
	if (!is.null(m <- x$markerdata)) {
		mis = attr(m, "missing")
		new_m = rbind(.prettyMarkers(m), matrix(mis, nrow=new_x$nInd-n, ncol=ncol(m)))
		new_x = setMarkers(new_x, new_m, missing=mis)
	}
	new_x
}


removeIndiv <- function(x, id) { #removes individual id and all his/her descendents. If the spouse is a founder, he/she will be removed as well.
	#founders without children after 'id' and 'desc' indivs are removed. The redundancy here does not matter.
	desc = descendents(x, id)
	leftover.spouses = setdiff(x$founders, as.numeric(x$pedigree[ -c(id, desc) , c('FID','MID')]))   
	
	new_df = as.data.frame(x)[-c(id, desc, leftover.spouses), ] 

	mis = ifelse(is.null(x$markerdata), 0, attr(x$markerdata, "missing"))
	linkdat(new_df, model=x$model, missing=mis)
}
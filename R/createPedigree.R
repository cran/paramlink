nuclearPed <- function(noffs, sex) {
    assert_that(.is.natural(noffs))
    if(missing(sex)) sex = rep.int(1, noffs)
    if(length(sex) < noffs) sex = rep(sex, length.out=noffs)
    p = cbind(ID=1:(2+noffs), FID=c(0, 0, rep.int(1, noffs)), MID=c(0, 0, rep.int(2, noffs)),
            SEX=c(1, 2, sex), AFF=1)
    linkdat(p, verbose=FALSE)
}


cousinsPed = function(degree, removal=0, degree2=NULL, child=FALSE) {
    # Creates a pedigree linking two n'th cousins k times removed, where n=degree, k=removal.
    # By default, removals are added on the right side, i.e. degree2 = degree + removal.
    # If degree2 is non-NULL, the removal parameter is ignored.
    assert_that(.is.natural0(degree), .is.natural0(removal), is.null(degree2) || .is.natural0(degree2))    
    if(is.null(degree2)) degree2 = degree + removal
    
    # Chain on the left side
    x = nuclearPed(1)
    for(i in seq_len(degree)) 
        x = addSon(x, x$nInd, verbose=F)
    
    # Chain on the right side
    x = addOffspring(x, father=1, mother=2, noffs=1, verbose=F)
    for(i in seq_len(degree2)) 
        x = addSon(x, x$nInd, verbose=F)
    x = swapSex(x, x$nInd, verbose=F)
    
    if(child) {
        cous = leaves(x)
        x = addOffspring(x, father=cous[1], mother=cous[2], noffs=1, verbose=F)
    }
    x
}


halfCousinsPed = function(degree, removal=0, degree2=NULL, child=FALSE) {
    # Creates a pedigree linking two n'th half cousins k times removed, where n=degree, k=removal.
    # By default, removals are added on the right side, i.e. degree2 = degree + removal.
    # If degree2 is non-NULL, the removal parameter is ignored.
    assert_that(.is.natural0(degree), .is.natural0(removal), is.null(degree2) || .is.natural0(degree2))
    if(is.null(degree2)) degree2 = degree + removal
    
    # Chain on the left side
    x = nuclearPed(1)
    for(i in seq_len(degree)) 
        x = addSon(x, x$nInd, verbose=F)
    
    # Chain on the right side
    x = addSon(x, 1, verbose=F)
    for(i in seq_len(degree2)) 
        x = addSon(x, x$nInd, verbose=F)
    x = swapSex(x, x$nInd, verbose=F)
    
    if(child) {
        cous = leaves(x)
        x = addOffspring(x, father=cous[1], mother=cous[2], noffs=1, verbose=F)
    }
    x
}

doubleCousins = function(degree1, degree2, removal1=0, removal2=0, child=FALSE) {
    # Creates a pedigree linking two individuals s.t their fathers are cousins of degree 'degree1',
    # and their mothers are cousins of degree 'degree2'.
    assert_that(.is.natural(degree1),.is.natural(degree2),.is.natural0(removal1),.is.natural0(removal2))
    
    # Pedigree linking the fathers
    x1 = cousinsPed(degree1-1, removal=removal1) 
    fathers = leaves(x1)
    x1 = swapSex(x1, fathers[2])
    
    # Pedigree linking the mothers 
    x2 = cousinsPed(degree2-1, removal=removal2)
    x2 = relabel(x2, new=1:x2$nInd + x1$nInd) # shifting labels
    mothers = leaves(x2)
    x2 = swapSex(x2, mothers[1])
    
    # Merging and adding children
    children = mothers[2] + 1:2
    p = rbind(as.matrix(x1, FALSE), as.matrix(x2, FALSE))
    p = rbind(p, c(1, children[1], fathers[1], mothers[1], 1,1), 
                 c(1, children[2], fathers[2], mothers[2], 2,1))
    x = linkdat(p, verbose=F)
    
    if(child)
        x = addOffspring(x, father=children[1], mother=children[2], noffs=1)
    x
}

# Wrapper for the most common case
doubleFirstCousins = function() doubleCousins(1,1)

quadHalfFirstCousins = function() {
    # Creates quad half fist cousins pedigree. NB: Does not draw well!
    p = matrix(c(
    1,0,0,1,1,
    2,0,0,2,1,
    3,0,0,1,1,
    4,0,0,2,1,
    5,1,2,1,1,
    6,3,4,2,1,
    7,1,4,1,1,
    8,3,2,2,1,
    9,5,6,1,1,
    10,7,8,2,1), byrow=T, nrow=10)
    linkdat(p, verbose=FALSE)
}

fullSibMating = function(generations) {
    # Creates a pedigree resulting from repeated brother-sister matings.
    assert_that(.is.natural(generations))
    x = nuclearPed(2, 1:2)
    for(i in seq_len(generations-1))
        x = addOffspring(x, father=x$nInd-1, mother=x$nInd, noffs=2, sex=1:2, verbose=F)
    x
}

halfSibStack = function(generations) {
    # Creates pedigree resulting from a breeding scheme where each generation adds two half brothers and a female founder. These become the parents of the half brothers in the next layer.
    assert_that(.is.natural(generations))
    x = linkdat(cbind(ID=1:5, FID=c(0,0,0,1,2), MID=c(0,0,0,3,3), SEX=c(1,1,2,1,1), AFF=1), verbose=FALSE)
    for(g in seq_len(generations)[-1]) {
        m = 3*g
        x = addOffspring(x, father=m-2, mother=m, noffs=1, verbose=F)
        x = addOffspring(x, father=m-1, mother=m, noffs=1, verbose=F)
    }
    x
}

#cousinCrossing = function(degrees, removals=0) {
#    # Creates pedigree resulting from a cousin-crossing breeding scheme.
#    # degrees: Vector of non-negative integers. (0 = siblings)
#    assert_that(is.numeric(degrees), length(degrees)>=1, all(degrees>=0), all(removals>=0))
#    removals = rep_len(removals, length(degrees))
#    x = cousinsPed(degrees[1], removals[1])
#    for(i in seq_along(degrees)[-1]) {
#        xn = x$nInd
#        xcousins = leaves(x)
#        y = cousinsPed(degrees[i], removals[i])
#        y = relabel(y, new = c(xcousins, seq.int(xn+1, length.out=y$nInd-2)))
#        x = mergePed(x, y, quick=TRUE)
#    }
#    x   
#}


#halfCousinStack = function(degrees, removals=0) {
#    # Creates pedigree by stacking layers of half cousins.
#    assert_that(is.numeric(degrees), length(degrees)>=1, all(degrees>=0), all(removals>=0))
#    removals = rep_len(removals, length(degrees))
#    x = halfCousinsPed(degrees[1], removals[1])
#    for(i in seq_along(degrees)[-1]) {
#        xn = max(x$orig.ids)
#        y = halfCousinsPed(degrees[i], removals[i])
#        y = relabel(y, new = xn+1:y$nInd)
#        xcousins = leaves(x)
#        x = mergePed(x, y, quick=TRUE)
#    }
#    x   
#}

mergePed <- function(x, y, quick=FALSE) {
    if(!is.null(x$markerdata) || !is.null(y$markerdata)) stop("Merging is only supported for pedigrees without marker data")
    ids = intersect(x$orig.ids, y$orig.ids)
    if(length(ids)==0) stop("Merging impossible: No common IDs")
    del = list(x = numeric(), y = numeric())
    for (i in ids) {
        if(.getSex(x,i) != .getSex(y,i)) stop(paste("Gender mismatch for individual", i))
        parx = parents(x, i)
        pary = parents(y, i)
        if(length(pary) == 0) del$y = c(del$y, i)
        else if(length(parx) == 0) del$x = c(del$x, i)
        else if(all(parx == pary)) del$y = c(del$y, i)
        else stop(paste("Parent mismatch for individual", i))
    }
    xx = as.matrix(x)[!x$orig.ids %in% del$x, ]
    yy = as.matrix(y)[!y$orig.ids %in% del$y, ]
    z = rbind(xx, yy)
    z[,'FAMID'] = x$famid # in case y$famid != x$famid
    
    if(quick) 
        return(restore_linkdat(z, attrs = attributes(xx), checkped=FALSE))
    
    # reorder to put parents above children (necessary when using IBDsim).
    N = nrow(z)
    i = 1
    while(i < N) {
        maxpar = max(match(z[i, c('FID','MID')], z[, 'ID'], nomatch=0))
        if(maxpar > i) {
            print("reordering!")
            z = z[c(seq_len(i-1), (i+1):maxpar, i, seq_len(N-maxpar)+maxpar), ]
        }
        else i = i+1
    }
    
    restore_linkdat(z, attrs = attributes(xx), checkped=TRUE)
}

######## DEPRECIATED ##########

# Replaced by cousins()
cousinPed <- function(degree) {
	stopifnot(degree>=0)
	if(degree==0) return(nuclearPed(noffs=2, sex=1:2))
	p = cbind(ID=1:4, FID=c(0,0,1,1), MID=c(0,0,2,2), SEX=c(1,2,1,1), AFF=1)
	for (n in 1:degree)
		p = rbind(p, c(4*n+1,0,0,2,1), c(4*n+2,0,0,2,1), c(4*n+3,4*n-1,4*n+1,1,1), c(4*n+4,4*n,4*n+2,1,1))
	p[nrow(p), 'SEX'] = 2
	linkdat(p, verbose=FALSE)
}

# Replaced by halfCousins()
halfCousinPed <- function(degree) {
	stopifnot(degree>=0)
	if(degree==0) p = cbind(ID=1:5, FID=c(0,0,0,1,1), MID=c(0,0,0,2,3), SEX=c(1,2,2,1,2), AFF=1)
	else {
		p = cbind(ID=1:3, FID=c(0,0,0), MID=c(0,0,0), SEX=c(1,2,2), AFF=1)
		for (n in seq_len(degree))
			p = rbind(p, c(4*n,4*n-4,4*n-3,1,1), c(4*n+1,0,0,2,1), c(4*n+2,4*n-2,4*n-1,1,1), c(4*n+3,0,0,2,1)) #add 1 generation: son in line 1, his wife, son in line 2, his wife
		dd = degree+1
		p = rbind(p, c(4*dd,4*dd-4,4*dd-3,1,1), c(4*dd+1,4*dd-2,4*dd-1,2,1)) #last generation - one boy, one girl.
		p[4,c(2,3)] = c(1,2); p[6,c(2,3)] = c(1,3) 
	}
	linkdat(p, verbose=FALSE)
}

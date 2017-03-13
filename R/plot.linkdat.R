
.as.kinship2_pedigree = function(x, deceased=numeric(0), aff2=NULL) {
    ped = relabel(x$pedigree, x$orig.ids) 
    AFF = ped[, 'AFF']; AFF[AFF == 0] = NA; AFF = AFF - 1 #kinship2:::plot.pedigree uses aff codes 1=aff, 0=unaff, NA=unknown
    if(!is.null(aff2)) {
        stopifnot(is.numeric(aff2), length(aff2)==x$nInd, all(aff2 %in% 0:2))
        aff2[aff2==0] = NA; aff2 = aff2 - 1
        AFF = cbind(aff1 = AFF, aff2 = aff2)
    }
    
    status = ifelse(x$orig.ids %in% deceased, 1, 0)
    pedigree(id=ped[,'ID'], dadid=ped[,'FID'], momid=ped[,'MID'], sex=ped[,'SEX'], affected=AFF, status=status)
}

plot.linkdat <- function(x, marker=NULL, alleles=NULL, sep="/", missing="-", 
                         skip.empty.genotypes=FALSE, id.labels=NULL, title=NULL, 
                         available=FALSE, col=1, deceased=numeric(0), starred=numeric(0), 
                         aff2=NULL, margins=c(0.6, 1, 4.1, 1), ...) {
    
    # Labels
    if(is.null(id.labels)) 
        if(!is.null(x$plot.labels)) {
            assert_that(length(x$plot.labels)==x$nInd)
            id.labels = x$plot.labels
        }
        else 
            id.labels = x$orig.ids
    else if(identical(id.labels, "num"))
        id.labels = x$orig.ids
    else if(!is.null(labs <- names(id.labels))) { 
        id.labels = labs[match(x$orig.ids, id.labels)]
        id.labels[is.na(id.labels)] = ""
    }
    
    if(!is.null(lb <- x$loop_breakers) && length(id.labels)==x$nInd) {
        origint = .internalID(x, lb[,1])
        copyint = .internalID(x, lb[,2])
        id.labels[copyint] = paste(id.labels[copyint], id.labels[origint], sep="=")
    }
    
    strid = rep(id.labels, length.out=x$nInd)
    
    # Add stars to labels
    starred = .internalID(x, starred)
    strid[starred] = paste0(strid[starred], "*")
    
    # Marker genotypes
    if (!is.null(marker)) {
        if (inherits(marker, "marker")) 
            m = list(marker)
        else if(is.list(marker) && all(sapply(marker, inherits, what="marker"))) 
            m = marker  #marker must be either a 'marker' object, a list of such, or an integer vector.
        else 
            m = x$markerdata[marker]
        gg = .prettyMarkers(m, alleles=alleles, sep=sep, missing=missing, singleCol=TRUE, sex=x$pedigree[, 'SEX'])
        geno = apply(gg, 1, paste, collapse="\n")
        if (skip.empty.genotypes) 
            geno[rowSums(do.call(cbind, m))==0] = ""
        if (is.null(strid) || !any(nzchar(strid))) 
            strid = geno 
        else 
            strid = paste(strid, geno, sep="\n")
    }
    oldmar = par(mar=margins)  # without this, par() does not equal 'margins'...(why??) Needed for centered title.
    
    # Colors
    cols = rep(col, length=x$nInd)
    
    # Coloring of available indivs
    avail_int = .internalID(x, x$available)
    if (is.logical(available) && isTRUE(available))
        cols[avail_int] = 2
    else if (tryCatch(length(col2rgb(available))==3, error=function(e) F)) # a single color
        cols[avail_int] = available
    
    # Special treatment for option "available=shaded"
    if (identical(available, "shaded")) {
        if(any(x$pedigree[,'AFF']==2)) 
            stop("The 'available=shaded' option cannot be used with disease pedigrees")
        if(any(c('angle','density') %in% names(list(...))))
            stop("Plot parameters 'angle' and 'density' cannot be used in combination with 'available=shaded'")
        x = swapAff(x, x$available)
        pedigree = .as.kinship2_pedigree(x, deceased=deceased, aff2=aff2)
        pdat = plot(pedigree, id=strid, col=cols, mar=margins, density=25, angle=45, ...)
    }
    else {
        pedigree = .as.kinship2_pedigree(x, deceased=deceased, aff2=aff2)
        pdat = plot(pedigree, id=strid, col=cols, mar=margins, ...)
    }
    
    # Add title
    if(!is.null(title)) 
        title(title)
    
    #par(oldmar)
    invisible(pdat)
}


plot.singleton = function(x, marker=NULL, alleles=NULL, sep="/", missing="-",
                          skip.empty.genotypes=FALSE, id.labels=NULL, title=NULL,
                          available=FALSE, col=1, deceased=numeric(0), starred=numeric(0), 
                          aff2=NULL, margins=c(8,0,0,0),  ...) {
    
    y = addParents(x, x$orig.ids, verbose=FALSE)
    p = plot.linkdat(y, marker=marker, alleles=alleles, sep=sep, missing=missing, 
                     skip.empty.genotypes=skip.empty.genotypes, id.labels=id.labels,
                     title=title, available=available, col=col, deceased=numeric(0), 
                     starred=starred, aff2=aff2, margins=c(margins[1],0,0,0), ...)
    usr = par('usr')
    rect(usr[1]-.1, p$y[3], usr[2]+.1, usr[4], border=NA, col="white")
    
    if(!is.null(title)) 
        title(title, line=-2.8)
}

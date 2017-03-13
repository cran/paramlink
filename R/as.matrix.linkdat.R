as.matrix.linkdat = function(x, include.attrs=TRUE, ...) {
    p = do.call(cbind, c(list(FAMID=x$famid, relabel(x$pedigree, x$orig.ids)), x$markerdata))
    if(include.attrs) {
        attr(p, "markerattr") = lapply(x$markerdata, attributes)
        attr(p, "available") = x$available
        attr(p, "model") = x$model
        attr(p, "plot.labels") = x$plot.labels
        attr(p, "orig.ids") = x$orig.ids # needed to restore plot.labels
    }
   p
}

restore_linkdat = function(x, attrs=NULL, checkped=TRUE) {
    if(is.null(attrs)) 
        attrs = attributes(x)
    y = linkdat(x[, 1:6, drop=F], model = attrs$model, checkped=checkped, verbose=FALSE)
    
    # Create marker objects
    markers = x[, -(1:6), drop=F]
    nMark = ncol(markers)/2
    if(nMark==0) markerdata_list = NULL
    else {
        markerattr = attrs$markerattr
        markerdata_list = lapply(seq_len(nMark), function(k) {
            m = markers[, c(2*k-1,2*k), drop=F]
            attributes(m) = c(markerattr[[k]][-1], list(dim=dim(m)))
            m
        })
        class(markerdata_list) = "markerdata"
    }
    
    y = setMarkers(y, markerdata_list)
    y = setAvailable(y, intersect(attrs$available, y$orig.ids))
    if(!is.null(pl <- attrs$plot.labels)) {
        y$plot.labels = pl[match(y$orig.ids, attrs$orig.ids)]
        y$plot.labels[is.na(y$plot.labels)] = ""
    }
    y
}

setPlotLabels = function(x, labels, ids=x$orig.ids) {
    assert_that(all(ids %in% x$orig.ids), length(labels)==length(ids))
    if(is.null(x$plot.labels)) 
        x$plot.labels = rep("", x$nInd)
    x$plot.labels[match(ids, x$orig.ids)] = labels
    x
}
    
.restore.linkdat = restore_linkdat
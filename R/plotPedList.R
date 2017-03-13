plotPedList = function(plot.arg.list, widths=NA, frames=T, frametitles=NULL, fmar=NA, newdev=F, dev.height=NA, dev.width=NA, ...) {
    
    plot.list.flattened = list()
    if(deduceFrames <- isTRUE(frames)) {frames = list(); k = 0}
    for(p in plot.arg.list) {
        if(is.linkdat(p)) p = list(p) # will now be included in next line
        if(is.linkdat.list(p))
            plot.list.flattened = c(plot.list.flattened, lapply(p, list))
        else {
            # if list of linkdat with plot arguments
            if(!is.linkdat(p[[1]])) {print(p); stop("First element is not a linkdat object.")}
            p = list(p)
            plot.list.flattened = append(plot.list.flattened, p)
        }
        if(deduceFrames) {
            group = (k+1):(k<-k+length(p))
            frames = append(frames, list(group))
        }
    }
    plot.arg.list = plot.list.flattened
    N = length(plot.arg.list)
    if(identical(widths, NA)) 
        widths = vapply(plot.arg.list, function(p) ifelse(is.singleton(p[[1]]), 1, 2.5), 1)
    else widths = rep_len(widths, N)
    maxGen = max(vapply(plot.arg.list, function(arglist) .generations(arglist[[1]]), 1))
    
    if(hasframetitles <- !is.null(frametitles)) 
        assert_that(is.character(frametitles), length(frametitles)==length(frames))
    
    extra.args = list(...)
    if(!'title' %in% names(extra.args))
        extra.args$title = ''
    
    defaultmargins = if(N>2) c(0,4,0,4) else c(0,2,0,2)
    
    plot.arg.list = lapply(plot.arg.list, function(arglist) {
        names(arglist)[1] = 'x'
        g = .generations(arglist$x)
        addMargin = 2*(maxGen-g+1)
        if(!'margins' %in% c(names(arglist), names(extra.args)))
            arglist$margins = defaultmargins + c(addMargin, 0, addMargin, 0)
        
        # additional arguments given in (...)
        for(parname in setdiff(names(extra.args), names(arglist))) 
            arglist[[parname]] = extra.args[[parname]]
        arglist
    })
    
    # title: this must be treated specially (in outer margins)
    titles = sapply(plot.arg.list, '[[', 'title')
    plot.arg.list = lapply(plot.arg.list, function(arglist) {arglist$title=""; arglist})
    hastitles = hasframetitles || any(titles != "")
    
    # frames: if list, check that each vector consists of consecutive integers, and no duplicates.
    if(is.list(frames)) {
        for(v in frames) 
            if(!identical(TRUE, all.equal(v, v[1]:v[length(v)])))
                stop(sprintf("Each element of 'frames' must consist of consecutive integers: %s", paste(v, collapse=',')))
        dup = anyDuplicated(unlist(frames))
        if(dup>0) stop(sprintf("Plot %d occurs twice in 'frames' list", dup))
    }
    
    
    # create layout of plot regions and plot!
    if(newdev) {
        if(is.na(dev.height)) dev.height = max(3, 1*maxGen) + .3*as.numeric(hastitles)
        if(is.na(dev.width)) dev.width=3*N
        dev.new(height=dev.height, width=dev.width, noRStudioGD=TRUE)
    }
    if(hastitles) 
        par(oma=c(0,0,3,0), xpd=NA)
    else 
        par(oma=c(0,0,0,0), xpd=NA)
    layout(rbind(1:N), widths=widths)
    for(arglist in plot.arg.list) 
        do.call(plot, arglist)
    
    # leftmost coordinate of each plot region (converted to value in [0,1]).
    ratios = c(0, cumsum(widths)/sum(widths))
    
    # add frames
    if (is.list(frames)) {
        midpoints = numeric()
        fstart_index = sapply(frames, function(v) v[1])
        fstop_index = sapply(frames, function(v) v[length(v)])
        ratio_start = ratios[fstart_index]
        ratio_stop = ratios[fstop_index + 1] # fordi 0 foerst
        midpoints = (ratio_start + ratio_stop)/2
        
        # margin (fmar): if NA, set to 5% of vertical height, but at most 0.25 inches.
        if(is.na(fmar)) fmar = min(0.05, 0.25/dev.size()[2])
        margPix = grconvertY(0, from='ndc', to="device")*fmar
        margXnorm = grconvertX(margPix, from='device', to="ndc")
        frame_start = grconvertX(ratio_start + margXnorm, from='ndc')
        frame_stop = grconvertX(ratio_stop - margXnorm, from='ndc')
        rect(xleft=frame_start, ybottom=grconvertY(1-fmar, from='ndc'), 
             xright=frame_stop, ytop=grconvertY(fmar, from='ndc'))
    }
    
    if(hasframetitles) for(i in 1:length(frames))
        mtext(frametitles[i], outer = TRUE, at=midpoints[i])
    else if(hastitles) for(i in 1:N)
        mtext(titles[i], outer = TRUE, at=(ratios[i]+ratios[i+1])/2)
    
}

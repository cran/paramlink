
exclusionPower = function(ped_claim, ped_true, ids, markerindex=NULL, alleles=NULL, afreq=NULL, known_genotypes=list(), Xchrom=FALSE, plot=TRUE) {
    if(is.linkdat(ped_claim))
        ped_claim = list(ped_claim)
    if(is.linkdat(ped_true))
        ped_true = list(ped_true)
   
    ids_claim = lapply(ped_claim, function(x) ids[ids %in% x$orig.ids])
    ids_true = lapply(ped_true, function(x) ids[ids %in% x$orig.ids])
   
    N_claim = length(ped_claim)
    N_true = length(ped_true)
    N = N_claim + N_true
   
    if(is.null(alleles)) {
        # Use markerdata of ped_claim and ped_true.
        # NB: No compatibility testing is done!! 
        partial_claim = lapply(ped_claim, function(p) p$markerdata[[markerindex]])
        partial_true = lapply(ped_true, function(p) p$markerdata[[markerindex]])
    }
    else {    
        if(length(alleles)==1) alleles=1:alleles
        if(Xchrom) chrom = 23 else chrom = NA
        partial_claim = lapply(1:N_claim, function(i) {
            x = ped_claim[[i]]
            m = marker(x, alleles=alleles, afreq=afreq, chrom=chrom)
            for(tup in known_genotypes) 
                if(tup[1] %in% x$orig.ids)
                    m = modifyMarker(x, m, ids=tup[1], genotype=tup[2:3])
            m})
        partial_true = lapply(1:N_true, function(i) {
            x = ped_true[[i]]
            m = marker(x, alleles=alleles, afreq=afreq, chrom=chrom)
            for(tup in known_genotypes) 
                if(tup[1] %in% x$orig.ids)
                    m = modifyMarker(x, m, ids=tup[1], genotype=tup[2:3])
            m})
    }
    
    if(isTRUE(plot) || plot=='plot_only') {
        op = par(oma=c(0,0,3,0), xpd=NA)
        widths = ifelse(sapply(c(ped_claim, ped_true), is.singleton), 1, 2)
        claim_ratio = sum(widths[1:N_claim])/sum(widths)
        layout(rbind(1:N), widths=widths)
        has_genotypes = length(known_genotypes)>0
        for(i in 1:N) {
            if (i<=N_claim) {x = ped_claim[[i]]; avail =  ids_claim[[i]]; mm = if(has_genotypes) partial_claim[[i]] else NULL}
            else {x = ped_true[[i - N_claim]]; avail = ids_true[[i - N_claim]]; mm = if(has_genotypes) partial_true[[i - N_claim]] else NULL}
            cols = ifelse(x$orig.ids %in% avail, 2, 1)
            plot(x, marker=mm, col=cols, margin=c(2,4,2,4), title="")
        }
        mtext("Claim", outer = TRUE, at=claim_ratio/2)
        mtext("True", outer = TRUE, at=.5+claim_ratio/2)
        rect(grconvertX(0.02, from='ndc'), grconvertY(0.02, from='ndc'),
            grconvertX(claim_ratio - 0.02, from='ndc'), grconvertY(.98, from='ndc'))
        rect(grconvertX(claim_ratio + 0.02, from='ndc'), grconvertY(0.02, from='ndc'),
            grconvertX(.98, from='ndc'), grconvertY(.98, from='ndc'))
        par(op)
        if(plot=='plot_only') return()
    }
   
    p.g = Reduce('%o%', lapply(which(lengths(ids_true)>0), function(i)
        oneMarkerDistribution(ped_true[[i]], ids=ids_true[[i]], partialmarker=partial_true[[i]], verbose=F)))
    I.g = Reduce('%o%', lapply(which(lengths(ids_claim)>0), function(i) 
        oneMarkerDistribution(ped_claim[[i]], ids=ids_claim[[i]], partialmarker=partial_claim[[i]], verbose=F)==0))
    sum(p.g*I.g)
}


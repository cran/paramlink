IBDest_getAlleleData = function(x, ids, markers=NULL) {
    stopifnot(is.linkdat.list(x), is.data.frame(ids))
    pednr1 = ids$pednr[1]
    pednr2 = ids$pednr[2]
    
    # Match IDs with orig.ids.
    int.id1 = match(ids$orig.ids[1], x[[pednr1]]$orig.ids)
    int.id2 = match(ids$orig.ids[2], x[[pednr2]]$orig.ids)
    
    # Collect alleles and frequencies in a matrix with 8 rows (first 4 alleles a,b,c,d followed by their freqs), one column per marker
    # Note that alleles are internal integers
    
    if(pednr1 == pednr2) {
        ped = x[[pednr1]]
        A = vapply(ped$markerdata[markers], function(m) {
            als = c(m[int.id1,], m[int.id2,])
            frq = rep_len(NA_real_, length(als))
            frq[als > 0] = attr(m, 'afreq')[als] # works, since 0's in als are ignored when indexing
            c(als, frq)
        }, FUN.VALUE=numeric(8))
    }
    else {
        ped1 = x[[pednr1]]
        ped2 = x[[pednr2]]
        A = vapply(markers, function(i) {
            m1 = ped1$markerdata[[i]]
            m2 = ped2$markerdata[[i]]
            als = c(m1[int.id1,], m2[int.id2,])
            frq = rep_len(NA_real_, length(als))
            frq[als > 0] = attr(m1, 'afreq')[als] # works, since 0's in als are ignored when indexing
            c(als, frq)
        }, FUN.VALUE=numeric(8))
    }
    return(A)
}


.IBDlikelihood = function(k,a,b,cc,d,pa,pb,pc,pd) {
    ### Vectorized function for computing kappa likelihoods, given genotypes for two related individuals
    # k: numeric of length 2 = (kappa0, kappa2)
    # a: vector of positive integers (allele 1 of individual 1)
    # b: vector of positive integers (allele 2 of individual 1)
    # cc: vector of positive integers (allele 1 of individual 2) (avoid overloading of 'c')
    # d: vector of positive integers (allele 2 of individual 2)
    # pa, pb, pc, pd: numeric vectors with frequencies of the above alleles.
    #
    # Note that all input vectors (except k) have the same length = #markers. 
    homoz1 = a==b
    homoz2 = cc==d
    mac = a==cc
    mbc = b==cc
    mad = a==d
    mbd = b==d
    g1.fr = 2^(!homoz1)*pa*pb
    g2.fr = 2^(!homoz2)*pc*pd 
    
    # Prob(g1, g2 | unrelated)
    UN = g1.fr * g2.fr
    
    # Prob(g1, g2 | parent-offspring)
    PO = (.5)^(homoz1+homoz2)*pa*pb*(pd*(mac+mbc) + pc*(mad+mbd))
    
    # Prob(g1, g2 | monozygotic twins)
    MZ = g1.fr * ((mac & mbd) | (mad & mbc))
    
    # return likelihoods (Thompson)
    k[1]*UN + (1-k[1]-k[2])*PO + k[2]*MZ
}

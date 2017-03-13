
pedigreeLoops = function(x) {
    dls = .descentPaths(x, 1:x$nInd, original.ids=FALSE)
    loops = list()
    for (id in 1:x$nInd) {
        if (length(dl <- dls[[id]])==1) next
        pairs = .comb2(length(dl))
        for (p in 1:nrow(pairs)) {
            p1 = dl[[pairs[p,1]]]; p2 = dl[[pairs[p,2]]]
            if (p1[2]==p2[2]) next
            inters = p1[match(p2,p1,0L)][-1]  #intersecting path members, excluding id
            if (length(inters)==0) next
            else {
                top=x$orig.ids[p1[1]]; bottom=x$orig.ids[inters[1]]
                pathA = p1[ seq_len(which(p1==inters[1]) - 2) + 1 ] #without top and bottom. Seq_len to deal with the 1:0 problem.
                pathB = p2[ seq_len(which(p2==inters[1]) - 2) + 1 ]
                loops = c(loops, list(list(top=top, bottom=bottom, pathA = x$orig.ids[pathA], pathB=x$orig.ids[pathB])))
            }
        }
    }
    unique(loops)
}

breakLoops = function(x, loop_breakers=NULL, verbose=TRUE) {
    if(is.singleton(x)) 
        stop("This function does not apply to singleton objects.")
    automatic = is.null(loop_breakers)
    if(automatic) {
        if(!x$hasLoops) return(x)
        loop_breakers = findLoopBreakers(x)
        if(length(loop_breakers)==0) {
            if(verbose) cat("Marriage loops detected, trying different selection method.\n")
            loop_breakers = findLoopBreakers2(x)
        }
    }
    
    if(length(loop_breakers)==0)
        stop("Loop breaking unsuccessful.")
    
    if(any(loop_breakers %in% x$orig.ids[x$founders])) 
        stop("Pedigree founders cannot be loop breakers.")

    if(verbose) 
        cat(ifelse(automatic,"Selected","User specified"), "loop breakers:", loop_breakers, '\n')
    
    pedm = as.matrix(x)#data.frame(x, famid=T, missing=0)
    attrs = attributes(pedm) #all attributes except 'dim'
    dup_pairs = x$loop_breakers #normally = NULL at this stage 
    for(id in loop_breakers) {
        dup_id = max(pedm[, 'ID']) + 1
        dup_pairs = rbind(dup_pairs, c(id, dup_id))
        intern = match(id, pedm[,'ID'])  #don't use .internalID here, since pedm changes all the time
        sex_col = pedm[intern, 'SEX'] + 2 # FID column if 'intern' is male; MID if female
        
        pedm = pedm[c(1:intern, intern, (intern+1):nrow(pedm)), ]
        pedm[intern+1, c('ID', 'FID','MID')] = c(dup_id, 0, 0)
        pedm[pedm[, sex_col] == id, sex_col] = dup_id
    }
    newx = restore_linkdat(pedm, attrs=attrs) 
    newx$loop_breakers = dup_pairs
    if(automatic) 
        return(breakLoops(newx, verbose=verbose))
    newx
}

findLoopBreakers = function(x) {
    loopdata = pedigreeLoops(x)
    # write each loop as vector exluding top/bottom
    loops = lapply(loopdata, function(lo) c(lo$pathA, lo$pathB))
    bestbreakers = numeric()
    while (length(loops) > 0) {
        # add the individual occuring in most loops
        best = which.max(tabulate(unlist(loops)))
        bestbreakers = c(bestbreakers, best)
        loops = loops[sapply(loops, function(vec) !best %in% vec)]
    }
    bestbreakers
}

findLoopBreakers2 = function(x) {
    if (!requireNamespace("igraph", quietly = TRUE))
        stop("This pedigree has marriage loops. The package 'igraph' must be installed for automatic selection of loop breakers.", call. = FALSE)

    ids = x$orig.ids
    N = max(ids)
    nonf = ids[x$nonfounders]
    breakers = numeric()
    
    ped2edge = function(p) {
        # input: ped-matrise UTEN founder-rader
        pp = cbind(p, paste(p[,'FID'], p[,'MID'], sep="+"))
        edge.children = pp[, c(4,1), drop=F]
        edge.marriage = unique.matrix(rbind(pp[, c(2,4), drop=F], pp[, c(3,4), drop=F]))
        rbind(edge.marriage, edge.children)
    }
    
    p = as.matrix(x, keep=F)[x$nonfounders, c('ID' ,'FID', 'MID'), drop=F]
    while(T) {
        g = igraph::graph_from_edgelist(ped2edge(p))
        loop = igraph::girth(g)$circle
        if(length(loop)==0) 
            break
        good.breakers = intersect(loop$name, nonf)
        if(length(good.breakers)==0)
            stop("This pedigree requires founders as loop breakers, which is not implemented in paramlink yet. Sorry!")
        b = good.breakers[1]
        breakers = c(breakers, b)
        
        b.is.parent = p == b & col(p) > 1
        p[b.is.parent] = as.character(N <- N+1)
    }
    breakers
}

tieLoops = function(x) {
    dups = x$loop_breakers
    if(is.null(dups) || nrow(dups)==0) {cat("No loops to tie\n"); return(x)}
    if(!all(dups %in% x$orig.ids)) stop("Something's wrong: Duplicated individuals no longer in pedigree.")
    pedm = as.matrix(x)
    attrs = attributes(pedm)

    origs = dups[,1]; copies = dups[,2]
    pedm = pedm[-match(copies, pedm[,'ID']), ]
    for(i in 1:length(origs)) {
        orig = origs[i]; copy = copies[i]
        sex = pedm[pedm[,'ID']==orig, 'SEX']
        pedm[pedm[, sex+2] == copy, sex+2] = orig 
    }
    restore_linkdat(pedm, attrs=attrs)
}

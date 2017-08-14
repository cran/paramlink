#' @export
all.equal.linkdat = function(target, current, ...) {
    res = TRUE
    if (!all.equal(class(target), class(current))) {
        cat("Class attributes are not equal\n")
        res = FALSE
    }
    if (target$nMark != current$nMark) {
        cat("Unequal numbers of markers:", target$nMark, "vs.", current$nMark, "\n")
        res = FALSE
    }
    names.t = names(target)
    names.c = names(current)
    if (!setequal(names.t, names.c)) {
        if (length(first_miss <- setdiff(names.c, names.t)) > 0) 
            cat("Missing slots in first object:", paste(first_miss, sep = ", "), "\n")
        if (length(sec_miss <- setdiff(names.t, names.c)) > 0) 
            cat("Missing slots in second object:", paste(sec_miss, sep = ", "), "\n")
        res = FALSE
    }
    
    # ID labels
    if (!setequal(target$orig.ids, current$orig.ids)) {
        cat("ID labels are not equal\n")
        res = FALSE
    }
    new_order = match(current$orig.ids, target$orig.ids)
    
    # Plot labels
    pl = current$plot.labels
    if (!is.null(pl) && !all(target$plot.labels[new_order] == pl)) {
        cat("Plot labels are not equal\n")
        res = FALSE
    }
    
    # Availability
    if (!setequal(target$available, current$available)) {
        cat("Unequal vectors of availability\n")
        res = FALSE
    }
    # Tree topologies
    ped_targ = relabel(target$pedigree, target$orig.ids)[new_order, , drop = F]
    ped_curr = relabel(current$pedigree, current$orig.ids)
    if (!identical(ped_curr, ped_targ)) {
        cat("Pedigree topologies are not equal\n")
        res = FALSE
    }
    
    if (!res) 
        return(res)
    
    if (target$nMark > 0) {
        mark_targ <- do.call(cbind, as.list(target$markerdata))[new_order, , drop = F]
        mark_curr <- do.call(cbind, as.list(current$markerdata))
        if (!isTRUE(all.equal(mark_targ, mark_curr))) {
            diffs = which(mark_targ != mark_curr, arr.ind = T)
            cat("Differences in the following markers:", sort(unique((diffs[, 2] + 1)%/%2)), 
                "\n")
            res = FALSE
        }
        markerattr_targ <- lapply(target$markerdata, attributes)
        markerattr_curr <- lapply(current$markerdata, attributes)
        if (!identical(markerattr_targ, markerattr_curr)) {
            diffattr = which(sapply(seq_along(markerattr_targ), function(i) !identical(markerattr_targ[[i]], 
                markerattr_curr[[i]])))
            cat("Difference in marker attributes for marker", diffattr, "\n")
            res = FALSE
        }
    }
    return(res)
}

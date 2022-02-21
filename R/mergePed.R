#' Merge two pedigrees
#'
#' This function merges two linkdat objects, joining them at the individuals
#' with equal ID labels. This is especially useful for building 'top-heavy'
#' pedigrees. Only linkdat objects without marker data are supported.
#'
#'
#' @param x,y \code{\link{linkdat}} objects
#' @param quick a single logical. If TRUE, no pedigree checks are performed, and
#'   the individual ordering may be unfortunate.
#' @return A \code{linkdat} object.
#'
#' @examples
#'
#' # Creating a trio where each parent have first cousin parents.
#' # (Alternatively, this could be built using many calls to addParents().)
#'
#' x = cousinPed(1)
#' x = addOffspring(x, father=7, mother=8, noffs=1, id=9)
#' x = addOffspring(x, father=9, mother=10, noffs=1, id=11)
#'
#' y = relabel(cousinPed(1), 101:108)
#' y = addOffspring(y, father=107, mother=108, noffs=1, sex=2, id=10)
#' y = addOffspring(y, father=9, mother=10, noffs=1, id=11)
#'
#' # Joining x and y at the common individuals 9,10,11:
#' z = mergePed(x,y)
#'
#' # plot all three pedigrees
#' op = par(mfrow = c(1,3))
#' plot(x); plot(y); plot(z)
#' par(op)
#'
#' @export
mergePed = function(x, y, quick = FALSE) {
    if (!is.null(x$markerdata) || !is.null(y$markerdata))
        stop("Merging is only supported for pedigrees without marker data")
    ids = intersect(x$orig.ids, y$orig.ids)
    if (length(ids) == 0)
        stop("Merging impossible: No common IDs")
    del = list(x = numeric(), y = numeric())
    for (i in ids) {
        if (.getSex(x, i) != .getSex(y, i))
            stop(paste("Gender mismatch for individual", i))
        parx = parents(x, i)
        pary = parents(y, i)
        if (length(pary) == 0)
            del$y = c(del$y, i) else if (length(parx) == 0)
            del$x = c(del$x, i) else if (all(parx == pary))
            del$y = c(del$y, i) else stop(paste("Parent mismatch for individual", i))
    }
    xx = as.matrix(x)[!x$orig.ids %in% del$x, ]
    yy = as.matrix(y)[!y$orig.ids %in% del$y, ]
    z = rbind(xx, yy)
    z[, "FAMID"] = x$famid  # in case y$famid != x$famid

    if (quick)
        return(restore_linkdat(z, attrs = attributes(xx), checkped = FALSE))

    # reorder to put parents above children (necessary when using IBDsim).
    z = .reorder_parents_before_children(z)
    restore_linkdat(z, attrs = attributes(xx), checkped = TRUE)
}


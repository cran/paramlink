#' Transfer marker data
#'
#' Transfer marker data between pedigrees (in the form of \code{\link{linkdat}}
#' objects). Both the source and target can be lists of linkdat and/or singleton
#' objects (these must have disjoint ID sets). Any previous marker data of the
#' target is overwritten.
#'
#' @param from a \code{\link{linkdat}} or \code{\link{singleton}} object, or a
#'   list of such objects.
#' @param to a \code{\link{linkdat}} or \code{\link{singleton}} object, or a
#'   list of such objects.
#' @return A \code{linkdat} object (or a list of such) similar to \code{to}, but
#' where all individuals also present in \code{from} have marker genotypes
#' copied over.  Any previous marker data is erased.
#'
#' @seealso \code{\link{linkdat}}
#'
#' @examples
#'
#' x = list(singleton(id=5), nuclearPed(noffs=2))
#' x = markerSim(x, N=5, alleles=1:5, verbose=FALSE, available=4:5)
#' y = nuclearPed(noffs=3)
#' y = transferMarkerdata(x, y)
#' stopifnot(all.equal(x[[1]], branch(y,5)))
#' stopifnot(all.equal(x[[2]], subset(y,1:4)))
#'
#' @export
transferMarkerdata = function(from, to) {
    assert_that(is.linkdat(from) || is.linkdat.list(from))
    assert_that(is.linkdat(to) || is.linkdat.list(to))

    if (is.linkdat(from) && is.linkdat(to))
        return(.transferMarkerdataSingle(from, to))
    if (is.linkdat(from) && is.linkdat.list(to))
        return(lapply(to, .transferMarkerdataSingle, from = from))
    if (is.linkdat.list(from) && is.linkdat(to)) {
        # start by transferring markers from the first in 'from'
        res = .transferMarkerdataSingle(from[[1]], to)
        b = as.matrix(res)

        # loop over the remaining and transfer
        for (from in from[-1]) {
            shared.ids = intersect(b[, "ID"], from$orig.ids)
            if (length(shared.ids) == 0)
                next
            a = as.matrix(from, include.attr = FALSE)
            b[match(shared.ids, b[, "ID"]), -(1:6)] = a[match(shared.ids, a[, "ID"]), -(1:6)]
        }
        y = restore_linkdat(b)
        available_int = rowSums(b[, -(1:6), drop = F]) > 0
        y = setAvailable(y, y$orig.ids[available_int])
        return(y)
    }
    if (is.linkdat.list(from) && is.linkdat.list(to))
        return(lapply(to, transferMarkerdata, from = from))
}


.transferMarkerdataSingle = function(from, to) {
    assert_that(is.linkdat(from), is.linkdat(to))
    if (from$nMark == 0) {
        warning("No markers to transfer.")
        return(to)
    }
    shared.ids = intersect(to$orig.ids, from$orig.ids)

    a = as.matrix(from)
    b = as.matrix(removeMarkers(to))
    allelematrix = a[, -(1:6), drop = FALSE]

    # create empty allele matrix, and copy rows of shared individuals
    allelematrix.new = matrix(0, ncol = ncol(allelematrix), nrow = to$nInd)
    allelematrix.new[match(shared.ids, b[, "ID"]), ] = allelematrix[match(shared.ids, a[, "ID"]),
        ]

    restore_linkdat(cbind(b, allelematrix.new), attrs = attributes(a))
}

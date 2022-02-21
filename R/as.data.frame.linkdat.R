#' linkdat to data.frame conversion
#'
#' Convert a linkdat object to data.frame for pretty printing.
#'
#' This function is mainly intended for pretty-printing \code{linkdat} objects
#' (for instance it is called by \code{print.linkdat}). For direct manipulation
#' of the pedigree and/or marker matrices, it is better to use
#' \code{\link{as.matrix.linkdat}}.
#'
#' @param x a \code{\link{linkdat}} object.
#' @param famid a logical indicating if the family identifier should be included
#'   as the first column.
#' @param markers a numeric indicating which markers should be included/printed.
#' @param alleles a character containing allele names, e.g.
#'   \code{alleles=c('A','B')}.
#' @param missing the character (of length 1) used for missing alleles. Defaults
#'   to '0'.
#' @param singleCol a logical: Should the two alleles for each marker be pasted
#'   into one column or kept in separate columns?
#' @param sep a single character to be used as allele separator if
#'   \code{singleCol=TRUE}.
#' @param \dots further arguments (not used).
#' @return A \code{data.frame}.
#' @seealso \code{\link{as.matrix.linkdat}}
#'
#' @examples
#'
#' x = linkdat(toyped)
#' x
#'
#' # Printing x as above is equivalent to:
#' as.data.frame(x, sep = '/', missing = '-', singleCol = TRUE)
#'
#' @export
as.data.frame.linkdat <- function(x, ..., famid = F, markers = seq_len(x$nMark), alleles = NULL,
    missing = NULL, singleCol = FALSE, sep = "") {
    p = relabel(x$pedigree, x$orig.ids)
    if (famid)
        p = cbind(FAMID = x$famid, p)

    if (!is.null(markers) && length(markers) > 0) {
        if (min(markers) < 1 || max(markers) > x$nMark)
            stop("Invalid marker number(s)")
        genotypes = .prettyMarkers(x$markerdata[markers], alleles = alleles, sep = sep, missing = missing,
            singleCol = singleCol, sex = p[, "SEX"])
    } else genotypes = matrix(numeric(0), nrow = x$nInd)

    data.frame(p, genotypes, stringsAsFactors = FALSE, ...)
}



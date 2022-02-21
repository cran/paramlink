#' Attach plot labels to a linkdat object
#'
#' This function attaches (or modifies) a character vector of plotting labels
#' for the pedigree members of a linkdat object. This is useful since only
#' numerical ID's are allowed in defining pedigrees in paramlink.
#'
#' @param x A linkdat object.
#' @param labels A character vector of the same length as \code{ids}.
#' @param ids A numeric vector of numerical IDs. Must be a subset of
#'   \code{x$orig.ids}.
#' @return A new linkdat object, differing from \code{x} only in
#'   \code{x$plot.labels}.
#' @seealso \code{\link{plot.linkdat}}
#'
#' @examples
#'
#' x = nuclearPed(1)
#' x = setPlotLabels(x, labels=c('Father', 'Mother', 'Son'))
#' plot(x)
#'
#' @export
setPlotLabels = function(x, labels, ids = x$orig.ids) {
    assert_that(all(ids %in% x$orig.ids), length(labels) == length(ids))
    if (is.null(x$plot.labels))
        x$plot.labels = rep("", x$nInd)
    x$plot.labels[match(ids, x$orig.ids)] = labels
    x
}

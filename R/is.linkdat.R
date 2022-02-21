#' Is an object a linkdat object?
#'
#' Functions for checking whether an object is a \code{\link{linkdat}} object, a
#' \code{\link{singleton}} or a list of such.
#'
#' Note that the \code{singleton} class inherits from \code{linkdat}, so if
#' \code{x} is a singleton, \code{is.linkdat(x)} returns TRUE.
#'
#' @param x Any \code{R} object.
#' @return For \code{is.linkdat}: TRUE if \code{x} is a linkdat (or singleton)
#'   object, and FALSE otherwise.\cr For \code{is.singleton}: TRUE if \code{x}
#'   is a singleton object, and FALSE otherwise.\cr For \code{is.linkdat.list}:
#'   TRUE if \code{x} is a list of linkdat/singleton objects.
#' @seealso \code{\link{linkdat}}
#' @examples
#'
#' x1 = nuclearPed(1)
#' x2 = singleton(1)
#' stopifnot(is.linkdat(x1), !is.singleton(x1),
#'           is.linkdat(x2), is.singleton(x2),
#'           is.linkdat.list(list(x1,x2)))
#'
#' @export
is.linkdat = function(x)
    inherits(x, "linkdat")

#' @rdname is.linkdat
#' @export
is.singleton = function(x)
    inherits(x, "singleton")

#' @rdname is.linkdat
#' @export
is.linkdat.list = function(x)
    isTRUE(is.list(x) && all(sapply(x, inherits, "linkdat")))

#' Relatedness coefficients
#'
#' Computes inbreeding coefficients for all pedigree members, and Jacquard's
#' condensed identity coefficients for any pair of members. These are simple
#' wrappers for functions in other packages or external programs.
#'
#' Both \code{inbreeding} and \code{kinship_coefs} are thin wrappers of
#' \code{\link[kinship2]{kinship}}. \code{jacquard2}, executes an external
#' call to the C program \code{IdCoefs} (Abney, 2009). For this to
#' function, \code{IdCoefs} must be installed on the computer (see link in the
#' References section below) and the executable placed in a folder included in
#' the PATH variable. The \code{jacquard2} wrapper works by writing the
#' necessary files to disk and calling \code{IdCoefs} via \code{\link{system}}.
#'
#' @param x a \code{\link{linkdat}} object.
#' @param ids a integer vector of length 2.
#' @param verbose a logical, indicating if messages from IdCoefs should be
#'   printed.
#' @param cleanup a logical: If TRUE, the pedfile and sample file created for
#'   the IdCoefs run are deleted automatically.
#' @return For \code{inbreeding}, a numerical vector with the inbreeding
#'   coefficients, with names according to the ID labels \code{x$orig.ids}.\cr
#'   For \code{kinship_coefs}, either a single numeric (if \code{ids} is a pair
#'   of pedigree members) or the whole kinship matrix, with \code{x$orig.ids} as
#'   dimnames.\cr For \code{jaquard} and \code{jaquard2}, a numerical vector of
#'   length 9 (in the standard order of Jacquard's condensed identity
#'   coefficients).
#' @seealso \code{\link[kinship2]{kinship}}
#' @references The \code{IdCoefs} program: Abney, Mark (2009). A graphical
#'   algorithm for fast computation of identity coefficients and generalized
#'   kinship coefficients. Bioinformatics, 25, 1561-1563.
#'   \url{http://home.uchicago.edu/~abney/abney_web/Software.html}
#'
#' @examples
#'
#' # Offspring of first cousins
#' x = cousinsPed(1, child=TRUE)
#' inb = inbreeding(x)
#' stopifnot(inb[9] == 1/16)
#'
#' # if ID labels are not 1:9, care must be taken in extracting correct elements.
#' set.seed(1357)
#' y = relabel(x, sample(1:9))
#' child = leaves(y)
#' inbreeding(y)[child] #wrong
#' inb = inbreeding(y)[as.character(child)] #correct
#' inb
#' # the inbreeding coeff of the child equals the kinship coeff of parents
#' kin = kinship_coefs(y, parents(y, child))
#' stopifnot(inb==kin, inb==1/16)
#'
#' @name relatednessCoeff
NULL

#' @rdname relatednessCoeff
#' @export
inbreeding = function(x) {
    ped = x$pedigree
    kin.matrix = kinship2::kinship(id = ped[, "ID"], dadid = ped[, "FID"], momid = ped[, "MID"])
    inb.coeff = numeric()
    inb.coeff[x$founders] = 0
    inb.coeff[x$nonfounders] = sapply(x$nonfounders, function(i) kin.matrix[ped[i, "FID"],
        ped[i, "MID"]])
    names(inb.coeff) = x$orig.ids
    inb.coeff
}

#' @rdname relatednessCoeff
#' @export
kinship_coefs = function(x, ids = NULL) {
    if (!is.null(ids))
        assert_that(length(ids) == 2, all(ids %in% x$orig.ids))
    ped = x$pedigree
    kin.matrix = kinship2::kinship(id = ped[, "ID"], dadid = ped[, "FID"], momid = ped[, "MID"])
    dimnames(kin.matrix) = list(x$orig.ids, x$orig.ids)
    if (is.null(ids))
        return(kin.matrix)
    kin.matrix[as.character(ids[1]), as.character(ids[2])]
}

#' @rdname relatednessCoeff
#' @export
jacquard = function(x, ids) {
    message("This function is no longer available. Use `ribd::identityCoefs()` instead")
    #if (!requireNamespace("identity", quietly = TRUE))
    #    stop("Package 'identity' must be install for this function to work.", call. = FALSE)
    #assert_that(length(ids) == 2, all(ids %in% x$orig.ids))
    #idsi = .internalID(x, ids)
    #ped = x$pedigree[, 1:3]
    #identity::identity.coefs(idsi, ped)[2, 3:11]
}


#' @rdname relatednessCoeff
#' @export
jacquard2 = function(x, ids, verbose = FALSE, cleanup = TRUE) {
    assert_that(length(ids) == 2, all(ids %in% x$orig.ids))
    x = .reorder_parents_before_children(x)
    ped = relabel(x$pedigree, x$orig.ids)[, 1:3] # relabel if orig.ids are 1..N (in some order)

    write.table(ped, file = "__paramlink2idcoefs__.ped", quote = F, row.names = F, col.names = F)
    write.table(ids, file = "__paramlink2idcoefs__.sample", quote = F, row.names = F, col.names = F)
    command = "idcoefs -p __paramlink2idcoefs__.ped -s __paramlink2idcoefs__.sample -o __paramlink2idcoefs__.output"
    run = suppressWarnings(system(command, intern = T))
    if (verbose)
        print(run)
    res = read.table("__paramlink2idcoefs__.output", as.is = T)
    if (cleanup)
        unlink(dir(pattern = "__paramlink2idcoefs__"))
    as.numeric(res[2, 3:11])
}

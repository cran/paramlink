#' Create simple pedigrees
#'
#' These are utility functions for creating some common pedigree structures as
#' \code{linkdat} objects.
#'
#' All individuals are created as unaffected. Use \code{\link{swapAff}} to edit
#' this (see Examples). Use \code{\link{swapSex}} to change gender of pedigree
#' members.
#'
#' The call \code{cousinsPed(degree=n, removal=k)} creates a pedigree with two
#' n'th cousins, k times removed. By default, removals are added on the right
#' side. To override this, the parameter \code{degree2} can be used to indicate
#' explicitly the number of generations on the right side of the pedigree. When
#' \code{degree2} is given \code{removal} is ignored. (Similarly for
#' \code{halfCousinsPed}.)
#'
#' The function \code{doubleCousins} creates two individuals whose fathers are
#' cousins (\code{degree1}, \code{removal1}) as well as their mothers
#' (\code{degree2}, \code{removal2}). For simplicity, a wrapper
#' \code{doubleFirstCousins} is provided for the most common case, double first
#' cousins. Finally \code{quadHalfFirstCousins} produces a pedigree with
#' quadruple half first cousins.
#'
#' \code{fullSibMating} crosses full sibs continuously for the indicated number
#' of generations.
#'
#' \code{halfSibStack} produces a breeding scheme where the two individuals in
#' the final generation are simultaneously half siblings and half n'th cousins,
#' where \code{n=1,...,generations}.
#'
#' \code{cousinPed} and \code{halfCousinPed} (written without the 's') are
#' depreciated functions kept for backwards compatibility. They create cousin
#' pedigrees, but without possibility for removals, and with a different
#' ordering than their replacements \code{cousinsPed} and \code{halfCousinsPed}.
#'
#'
#' @param noffs A positive integer, the number of offspring in the nuclear
#'   family.
#' @param sex A vector of length \code{noffs}; indicating the genders (1=male,
#'   2=female) of the offspring. If missing, all offspring are taken to be
#'   males.
#' @param degree,degree1,degree2 Non-negative integers, indicating the degree of
#'   cousin-like relationships: 0=siblings, 1=first cousins; 2=second cousins,
#'   a.s.o. See Details and Examples.
#' @param removal,removal1,removal2 Non-negative integers, indicating removals
#'   of cousin-like relationships. See Details and Examples.
#' @param child A logical: Should an inbred child be added to the two cousins?
#' @param generations A positive integer indicating the number of crossings.
#' @return A \code{\link{linkdat}} object.
#' @seealso \code{\link{swapAff}}, \code{\link{swapSex}},
#'   \code{\link{removeIndividuals}}, \code{\link{addOffspring}},
#'   \code{\link{relabel}}
#'
#' @examples
#'
#' # A nuclear family with 2 boys and 3 girls,
#' # where the father and the two boys are affected.
#' x = nuclearPed(noffs=5, sex=c(1,1,2,2,2))
#' x = swapAff(x, ids=c(1,3,4))
#'
#' # Half sibs:
#' halfCousinsPed(degree=0)
#'
#' # Grand aunt:
#' cousinsPed(degree=0, removal=2)
#'
#' # Second cousins once removed.
#' cousinsPed(degree=2, removal=1)
#'
#' # Again second cousins once removed,
#' # but with the 'removal' on the left side.
#' cousinsPed(degree=3, degree2=2)
#'
#' # A child of first cousin parents.
#' cousinsPed(degree=1, child=TRUE)
#'
#' # Consecutive brother-sister matings.
#' fullSibMating(3)
#'
#' # Simultaneous half siblings and half first cousins
#' halfSibStack(2)
#'
#' # Double first cousins
#' doubleFirstCousins()
#'
#' # Quadruple half first cousins
#' # Weird plotting behaviour for this pedigree.
#' x = quadHalfFirstCousins()
#' #plot(x)
#'
#' @name pedCreate
NULL

#' @rdname pedCreate
#' @export
nuclearPed = function(noffs, sex) {
    assert_that(.is.natural(noffs))
    if (missing(sex))
        sex = rep.int(1, noffs) else assert_that(length(sex) <= noffs)
    if (length(sex) < noffs)
        sex = rep(sex, length.out = noffs)

    p = cbind(ID = 1:(2 + noffs), FID = c(0, 0, rep.int(1, noffs)), MID = c(0, 0, rep.int(2,
        noffs)), SEX = c(1, 2, sex), AFF = 1)
    linkdat(p, verbose = FALSE)
}

#' @rdname pedCreate
#' @export
cousinsPed = function(degree, removal = 0, degree2 = NULL, child = FALSE) {
    # Creates a pedigree linking two n'th cousins k times removed, where n=degree, k=removal.
    # By default, removals are added on the right side, i.e. degree2 = degree + removal.  If
    # degree2 is non-NULL, the removal parameter is ignored.
    assert_that(.is.natural0(degree), .is.natural0(removal), is.null(degree2) || .is.natural0(degree2))
    if (is.null(degree2))
        degree2 = degree + removal

    # Chain on the left side
    x = nuclearPed(1)
    for (i in seq_len(degree)) x = addSon(x, x$nInd, verbose = F)

    # Chain on the right side
    x = addOffspring(x, father = 1, mother = 2, noffs = 1, verbose = F)
    for (i in seq_len(degree2)) x = addSon(x, x$nInd, verbose = F)
    x = swapSex(x, x$nInd, verbose = F)

    if (child) {
        cous = leaves(x)
        x = addOffspring(x, father = cous[1], mother = cous[2], noffs = 1, verbose = F)
    }
    x
}

#' @rdname pedCreate
#' @export
halfCousinsPed = function(degree, removal = 0, degree2 = NULL, child = FALSE) {
    # Creates a pedigree linking two n'th half cousins k times removed, where n=degree,
    # k=removal.  By default, removals are added on the right side, i.e. degree2 = degree +
    # removal.  If degree2 is non-NULL, the removal parameter is ignored.
    assert_that(.is.natural0(degree), .is.natural0(removal), is.null(degree2) || .is.natural0(degree2))
    if (is.null(degree2))
        degree2 = degree + removal

    # Chain on the left side
    x = nuclearPed(1)
    for (i in seq_len(degree)) x = addSon(x, x$nInd, verbose = F)

    # Chain on the right side
    x = addSon(x, 1, verbose = F)
    for (i in seq_len(degree2)) x = addSon(x, x$nInd, verbose = F)
    x = swapSex(x, x$nInd, verbose = F)

    if (child) {
        cous = leaves(x)
        x = addOffspring(x, father = cous[1], mother = cous[2], noffs = 1, verbose = F)
    }
    x
}

#' @rdname pedCreate
#' @export
doubleCousins = function(degree1, degree2, removal1 = 0, removal2 = 0, child = FALSE) {
    # Creates a pedigree linking two individuals s.t their fathers are cousins of degree
    # 'degree1', and their mothers are cousins of degree 'degree2'.
    assert_that(.is.natural(degree1), .is.natural(degree2), .is.natural0(removal1), .is.natural0(removal2))

    # Pedigree linking the fathers
    x1 = cousinsPed(degree1 - 1, removal = removal1)
    fathers = leaves(x1)
    x1 = swapSex(x1, fathers[2])

    # Pedigree linking the mothers
    x2 = cousinsPed(degree2 - 1, removal = removal2)
    x2 = relabel(x2, new = 1:x2$nInd + x1$nInd)  # shifting labels
    mothers = leaves(x2)
    x2 = swapSex(x2, mothers[1])

    # Merging and adding children
    children = mothers[2] + 1:2
    p = rbind(as.matrix(x1, FALSE), as.matrix(x2, FALSE))
    p = rbind(p, c(1, children[1], fathers[1], mothers[1], 1, 1), c(1, children[2], fathers[2],
        mothers[2], 2, 1))
    x = linkdat(p, verbose = F)

    if (child)
        x = addOffspring(x, father = children[1], mother = children[2], noffs = 1)
    x
}

#' @rdname pedCreate
#' @export
doubleFirstCousins = function()
    # Wrapper for the most common case of doubleCousins()
    doubleCousins(1, 1)

#' @rdname pedCreate
#' @export
quadHalfFirstCousins = function() {
    # Creates quad half fist cousins pedigree. NB: Does not draw well!
    p = matrix(c(1, 0, 0, 1, 1, 2, 0, 0, 2, 1, 3, 0, 0, 1, 1, 4, 0, 0, 2, 1, 5, 1, 2, 1, 1,
        6, 3, 4, 2, 1, 7, 1, 4, 1, 1, 8, 3, 2, 2, 1, 9, 5, 6, 1, 1, 10, 7, 8, 2, 1), byrow = T,
        nrow = 10)
    linkdat(p, verbose = FALSE)
}

#' @rdname pedCreate
#' @export
fullSibMating = function(generations) {
    # Creates a pedigree resulting from repeated brother-sister matings.
    assert_that(.is.natural(generations))
    x = nuclearPed(2, 1:2)
    for (i in seq_len(generations - 1)) x = addOffspring(x, father = x$nInd - 1, mother = x$nInd,
        noffs = 2, sex = 1:2, verbose = F)
    x
}

#' @rdname pedCreate
#' @export
halfSibStack = function(generations) {
    # Creates pedigree resulting from a breeding scheme where each generation adds two half
    # brothers and a female founder.  These become the parents of the half brothers in the next
    # layer.
    assert_that(.is.natural(generations))
    x = linkdat(cbind(ID = 1:5, FID = c(0, 0, 0, 1, 2), MID = c(0, 0, 0, 3, 3), SEX = c(1,
        1, 2, 1, 1), AFF = 1), verbose = FALSE)
    for (g in seq_len(generations)[-1]) {
        m = 3 * g
        x = addOffspring(x, father = m - 2, mother = m, noffs = 1, verbose = F)
        x = addOffspring(x, father = m - 1, mother = m, noffs = 1, verbose = F)
    }
    x
}




######## DEPRECIATED ##########

#' @rdname pedCreate
#' @export
cousinPed <- function(degree) {
    # Replaced by cousins()
    stopifnot(degree >= 0)
    if (degree == 0)
        return(nuclearPed(noffs = 2, sex = 1:2))
    p = cbind(ID = 1:4, FID = c(0, 0, 1, 1), MID = c(0, 0, 2, 2), SEX = c(1, 2, 1, 1), AFF = 1)
    for (n in 1:degree) p = rbind(p, c(4 * n + 1, 0, 0, 2, 1), c(4 * n + 2, 0, 0, 2, 1), c(4 *
        n + 3, 4 * n - 1, 4 * n + 1, 1, 1), c(4 * n + 4, 4 * n, 4 * n + 2, 1, 1))
    p[nrow(p), "SEX"] = 2
    linkdat(p, verbose = FALSE)
}

#' @rdname pedCreate
#' @export
halfCousinPed <- function(degree) {
    # Replaced by halfCousins()
    stopifnot(degree >= 0)
    if (degree == 0)
        p = cbind(ID = 1:5, FID = c(0, 0, 0, 1, 1), MID = c(0, 0, 0, 2, 3), SEX = c(1, 2, 2,
            1, 2), AFF = 1) else {
        p = cbind(ID = 1:3, FID = c(0, 0, 0), MID = c(0, 0, 0), SEX = c(1, 2, 2), AFF = 1)
        for (n in seq_len(degree)) p = rbind(p, c(4 * n, 4 * n - 4, 4 * n - 3, 1, 1), c(4 *
            n + 1, 0, 0, 2, 1), c(4 * n + 2, 4 * n - 2, 4 * n - 1, 1, 1), c(4 * n + 3, 0, 0,
            2, 1))  #add 1 generation: son in line 1, his wife, son in line 2, his wife
        dd = degree + 1
        p = rbind(p, c(4 * dd, 4 * dd - 4, 4 * dd - 3, 1, 1), c(4 * dd + 1, 4 * dd - 2, 4 *
            dd - 1, 2, 1))  #last generation - one boy, one girl.
        p[4, c(2, 3)] = c(1, 2)
        p[6, c(2, 3)] = c(1, 3)
    }
    linkdat(p, verbose = FALSE)
}

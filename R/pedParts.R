#' Pedigree subsets
#'
#' Utility functions for 'linkdat' objects, mainly for extracting various
#' pedigree information.
#'
#' @param x a \code{\link{linkdat}} object. In \code{related.pairs} possibly a
#'   list of \code{linkdat} objects.
#' @param id a numerical ID label.
#' @param original.id a logical indicating whether 'id' refers to the original
#'   ID label or the internal labeling.
#' @param degree a non-negative integer.
#' @param removal a non-negative integer
#' @param half a logical or NA. If TRUE (resp FALSE), only half (resp. full)
#'   siblings/cousins/nephews/nieces are returned. If NA, both categories are
#'   included.
#' @param relation one of the words (possibly truncated) \code{parents},
#'   \code{siblings}, \code{grandparents}, \code{nephews_nieces},
#'   \code{cousins}, \code{spouses}, \code{unrelated}.
#' @param interfam one of the words (possibly truncated) \code{none},
#'   \code{founders} or \code{all}, specifying which interfamiliar pairs should
#'   be included as unrelated in the case where \code{x} is a list of several
#'   pedigrees. If \code{none}, only intrafamiliar pairs are considered; if
#'   \code{founders} all interfamiliar pairs of (available) founders are
#'   included; if \code{all}, all interfamiliar (available) pairs are included.
#' @param available a logical, if TRUE only pairs of available individuals are
#'   returned.
#' @param ... further parameters
#'
#' @return For \code{ancestors(x,id)}, a vector containing the ID's of all
#'   ancestors of the individual \code{id}.  For \code{descendants(x,id)}, a
#'   vector containing the ID's of all descendants (i.e. children,
#'   grandchildren, a.s.o.) of individual \code{id}.
#'
#'   The functions \code{cousins}, \code{grandparents}, \code{nephews_nieces},
#'   \code{offspring}, \code{parents}, \code{siblings}, \code{spouses},
#'   \code{unrelated}, each returns an integer vector containing the ID's of all
#'   pedigree members having the specified relationship with \code{id}.
#'
#'   For \code{related.pairs} a matrix with two columns. Each row gives of a
#'   pair of pedigree members with the specified relation. If the input is a
#'   list of multiple pedigrees, the matrix entries are characters of the form
#'   'X-Y' where X is the family ID and Y the individual ID of the person.
#'
#'   For \code{leaves}, a vector of IDs containing all pedigree members without
#'   children.
#'
#' @examples
#'
#' p = cbind(ID=2:9, FID=c(0,0,2,0,4,4,0,2), MID=c(0,0,3,0,5,5,0,8),
#'         SEX=c(1,2,1,2,1,2,2,2), AFF=c(2,1,2,1,2,1,1,2))
#' x = linkdat(p)
#' stopifnot(setequal(spouses(x, 2), c(3,8)),
#'           setequal(offspring(x, 2), c(4,9)),
#'           setequal(descendants(x, 2), c(4,6,7,9)),
#'           setequal(leaves(x), c(6,7,9)))
#'
#' # Creating a loop and detecting it with 'pedigreeLoops'
#' # (note that we get two loops, one for each inbred child):
#' loopx = addOffspring(x, father=4, mother=9, noffs=2)
#' lps = pedigreeLoops(loopx)
#' stopifnot(lps[[1]]$top == 2, setequal(sapply(lps, '[[', 'bottom'), 10:11))
#'
#' # We add genotypes for a single SNP marker and compute a LOD score under a dominant model.
#' loopx = setMarkers(loopx, cbind(1,c(2,1,2,1,2,1,1,2,1,1)))
#' loopx = setModel(loopx, 1)
#'
#' # Loops are automatically broken in lod():
#' LOD1 = lod(loopx, theta=0.1)
#' stopifnot(round(LOD1, 3) == 1.746)
#'
#' # Or we can break the loop manually before computing the LOD:
#' loopfree = breakLoops(loopx, loop_breaker=4)
#' LOD2 = lod(loopfree, theta=0.1)
#' stopifnot(all.equal(loopx, tieLoops(loopfree)))
#' stopifnot(all.equal(LOD1, LOD2))
#'
#' @name pedParts
NULL


#' @rdname pedParts
#' @export
offspring = function(x, id, original.id = TRUE) {
    if (original.id)
        id = .internalID(x, id)
    p = x$pedigree
    offs_rows = p[, 1 + p[id, "SEX"]] == id
    if (original.id)
        x$orig.ids[offs_rows] else (1:x$nInd)[offs_rows]
}

#' @rdname pedParts
#' @export
spouses = function(x, id, original.id = TRUE) {
    # Returns a vector containing all individuals sharing offspring with <id>.
    internal_id = ifelse(original.id, .internalID(x, id), id)
    p = x$pedigree
    offs_rows = p[, 1 + p[internal_id, "SEX"]] == internal_id
    spou = unique.default(p[offs_rows, 4 - p[internal_id, "SEX"]])  # sex=1 -> column 3; sex=2 -> column 2.
    if (original.id)
        return(x$orig.ids[spou]) else return(spou)
}

#' @rdname pedParts
#' @export
related.pairs = function(x, relation = c("parents", "siblings", "grandparents", "nephews_nieces",
    "cousins", "spouses", "unrelated"), available = F, interfam = c("none",
    "founders", "all"), ...) {
    relation = match.arg(relation)
    interfam = match.arg(interfam)
    func = function(...) get(relation)(...)

    if (is.linkdat.list(x)) {
        res = do.call(rbind, lapply(x, function(xx)
            related.pairs(xx, relation, available, ...)
        ))
        if (relation == "unrelated" && interfam != "none") {
            avail = lapply(x, function(xx) {
                ids = if (available)
                  xx$available else xx$orig.ids
                if (interfam == "founders")
                  ids = intersect(ids, xx$orig.ids[xx$founders])
                if (length(ids) == 0)
                  return(NULL)
                ids
            })
            avail = avail[!sapply(avail, is.null)]
            fampairs = data.frame(t(.comb2(length(avail))))  # enable use of lapply below
            interfam = do.call(rbind, lapply(fampairs, function(p) fast.grid(avail[p])))
            res = rbind(res, interfam)
        }
        return(res)
    }

    res = NULL
    for (i in 1:x$nInd) {
        rels = func(x, i, original.id = F, ...)
        rels = rels[rels != i]
        res = rbind(res, cbind(rep.int(i, length(rels)), rels, deparse.level = 0))
    }
    res[res[, 1] > res[, 2], ] = res[res[, 1] > res[, 2], 2:1]
    res = unique(res)
    if (available) {
        avail = .internalID(x, x$available)
        res = res[res[, 1] %in% avail & res[, 2] %in% avail, , drop = F]
    }

    # return matrix with original IDs
    res[] = x$orig.ids[res]
    res
}

#' @rdname pedParts
#' @export
unrelated = function(x, id, original.id = TRUE) {
    if (!original.id)
        id = x$orig.ids[id]
    ancs = c(id, ancestors(x, id))
    rel = unique.default(unlist(lapply(ancs, function(a) c(a, descendants(x, a, original.id = TRUE)))))
    unrel = setdiff(x$orig.ids, rel)
    if (!original.id)
        unrel = .internalID(x, unrel)
    unrel
}


#' @rdname pedParts
#' @export
leaves = function(x) {
    p = as.matrix(x, FALSE)
    .mysetdiff(p[, "ID", drop = F], p[, c("FID", "MID")])
}

#' @rdname pedParts
#' @export
parents = function(x, id, original.id = TRUE) {
    grandparents(x, id, degree = 1, original.id = original.id)
}

#' @rdname pedParts
#' @export
grandparents = function(x, id, degree = 2, original.id = TRUE) {
    if (original.id)
        id = .internalID(x, id)
    p = x$pedigree
    gp = id
    for (i in seq_len(degree)) gp = p[gp, 2:3]
    if (original.id)
        x$orig.ids[gp] else (1:x$nInd)[gp]
}

#' @rdname pedParts
#' @export
siblings = function(x, id, half = NA, original.id = TRUE) {
    if (original.id)
        id = .internalID(x, id)
    p = x$pedigree
    fa = p[id, "FID"]
    mo = p[id, "MID"]
    if (fa == 0 && mo == 0)
        return(numeric())
    samefather = p[, "FID"] == fa
    samemother = p[, "MID"] == mo
    sib_rows = if (is.na(half))
        samefather | samemother else if (half)
        xor(samefather, samemother) else samefather & samemother
    sib_rows[id] = FALSE
    if (original.id)
        x$orig.ids[sib_rows] else (1:x$nInd)[sib_rows]
}

#' @rdname pedParts
#' @export
cousins = function(x, id, degree = 1, removal = 0, half = NA, original.id = TRUE) {
    if (original.id)
        id = .internalID(x, id)
    gp = grandparents(x, id, degree = degree, original.id = FALSE)
    uncles = unique.default(unlist(lapply(gp, siblings, x = x, half = half, original.id = FALSE)))
    cous = uncles
    for (i in seq_len(degree + removal)) cous = unique.default(unlist(lapply(cous, offspring,
        x = x, original.id = FALSE)))
    if (original.id)
        cous = x$orig.ids[cous]
    cous
}

#' @rdname pedParts
#' @export
nephews_nieces = function(x, id, removal = 1, half = NA, original.id = TRUE) {
    cousins(x, id, degree = 0, removal = removal, half = half, original.id = original.id)
}

#' @rdname pedParts
#' @export
ancestors = function(x, id) {
    # climbs up the pedigree storing parents iteratively. (Not documented: Accepts id of length
    # > 1)
    if (is.linkdat(x)) {
        p = x$pedigree
        orig_ids = x$orig.ids
        ids_int = .internalID(x, id)
    } else if (is.matrix(x) && isTRUE(all(c("ID", "FID", "MID") %in% colnames(x)))) {
        p = x
        orig_ids = p[, "ID"]
        ids_int = match(id, orig_ids)
    } else stop("x must be either a linkdat object or a matrix whose colnames include 'ID', 'FID' and 'MID'")
    p = relabel(p, 1:nrow(p))

    ancest = numeric(0)
    up1 = as.numeric(p[ids_int, c("FID", "MID")])
    up1 = up1[up1 > 0 & up1 <= nrow(p)]  #NB: Avoids pedigree errors without warning! Should be caught in .checkped anyway
    up1 = up1[!duplicated.default(up1)]
    while (length(up1) > 0) {
        ancest = c(ancest, up1)
        up1 = .mysetdiff(as.numeric(p[up1, c("FID", "MID")]), ancest)
    }
    ancest = sort.int(ancest[(ancest != 0) & !duplicated(ancest)])
    return(orig_ids[ancest])
}

#' @rdname pedParts
#' @export
descendants = function(x, id, original.id = TRUE) {
    if(original.id)
        id = .internalID(x, id)
    ped = x$pedigree
    ID = ped[, 'ID']
    F = ped[, 'FID']
    M = ped[, 'MID']

    desc = numeric()
    nextoffs = id
    while(length(nextoffs)) {
        nextoffs = ID[F %in% nextoffs | M %in% nextoffs]
        desc = c(desc, nextoffs)
    }
    desc = sort.int(unique.default(desc))
    if (original.id)
        desc = x$orig.ids[desc]
    desc
}

#' Pairwise common ancestors
#'
#' Computes a matrix A whose entry A[i,j] is TRUE if pedigree members i and j have a common ancestor, and FALSE otherwise.
#'
#' @param x a \code{\link{linkdat}} object.
#'
#' @examples
#'
#' x = fullSibMating(3)
#' A = hasCA(x)
#' stopifnot(A[1,1], !A[1,2], all(A[3:8, 3:8]))
#'
#' @export
hasCA = function(x) {
    if(!all(x$orig.ids - 1:x$nInd == 0))
        stop("This is currently only implemented for pedigrees with ordering 1,2,...")
    A = matrix(F, ncol=x$nInd, nrow=x$nInd)
    for(i in x$founders) {
        # vector of all descendants of i, including i
        desc = c(i, descendants(x,i))
        A[.my.grid(rep(list(desc), 2))] = T
    }
    A
}

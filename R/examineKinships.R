#' Check pedigree for relationship errors
#' 
#' This function provides a convenient way to check for pedigree errors in a 
#' linkage project or other situations where marker data is available for several 
#' members. The function estimates IBD coefficients for all indicated
#' pairs in the pedigree(s) and produces a color-coded plot where wrong
#' relationships are easy to spot.
#' 
#' This function works with SNP markers only (for now).
#' 
#' @param x A \code{\link{linkdat}} object, or a list of such.
#' @param who A character vector of one or more of the words 'parents',
#'   'siblings', 'grandparents', 'cousins', 'distant' and 'unrelated'. Two short
#'  forms are possible: 'all' (all of the above) and 'close' (all of the above
#'  except 'distant' and 'unrelated'.)
#' @param interfam A character; either 'founders', 'none' or 'all', indicating
#'  which interfamiliar pairs of individuals should be included. Only relevant
#'  if \code{x} is a list of several \code{linkdat} objects.
#' @param makeplot A logical.
#' @param showRelationships A character vector passed on to the
#'  \code{relationships} argument of \code{IBDtriangle}.
#' @param \dots Other arguments passed on to \code{IBDestimate}.
#'
#' @return A list of matrices (one for each relation category) with IBD estimates.
#' @author Magnus Dehli Vigeland
#' @seealso \code{\link{IBDestimate}}, \code{\link{IBDtriangle}}
#'
#' @examples
#' 
#' x = cousinsPed(1)
#' x = simpleSim(x, 500, alleles=1:2)
#' examineKinships(x)
#'
#' # Pretend we didn't know the borthers were related
#' x1 = branch(x, 3)
#' x2 = branch(x, 6)
#' x2$famid = 2
#'
#' # Notice the grey dot close to the 'S' (siblings)
#' examineKinships(list(x1, x2))
#'
#' @export
examineKinships = function(x, who = "all", interfam = c("founders", "none", "all"), makeplot = T, 
    showRelationships = c("UN", "PO", "MZ", "S", "H,U,G", "FC", "SC"), ...) {
    
    linkdatList = is.linkdat.list(x)
    if (linkdatList && any(duplicated(sapply(x, function(xx) xx$famid)))) 
        stop("When input is a list of linkdat objects, they must have different family ID (famid)")
    
    rel.groups = c("parents", "siblings", "grandparents", "cousins", "distant", "unrelated")
    if (identical(who, "all")) {
        groups = rel.groups
        if (linkdatList) 
            groups = groups[groups != "distant"]
    } else if (identical(who, "close")) 
        groups = rel.groups[1:4] else groups = sapply(who, match.arg, rel.groups[-5])
    
    
    if (makeplot) {
        # oldpar = par(no.readonly = TRUE)
        if (is.linkdat(x)) {
            g = .generations(x)
            marg = max(1, 12 - 2 * g)
            layout(rbind(1:2), widths = c(0.4, 0.6))
            plot(x, av = T, cex = 0.8, title = "", margins = c(marg, 1, marg, 1))
        }
        IBDtriangle(relationships = showRelationships)
        leg_txt = c("Parent-offspring", "Siblings", "Halfsib/Uncle/Grand", "1st cousins", "Distant", 
            "Unrelated")[rel.groups %in% groups]
        leg_col = c(2, 4, 3, 6, "cyan", "darkgray")[rel.groups %in% groups]
        legend("topright", title = " According to pedigree:", title.adj = 0, legend = leg_txt, 
            col = leg_col, pch = 16)
    }
    
    dat = .ibdPrep(x)
    p = s = g = co = d = u = NULL
    
    if ("parents" %in% groups) {
        P = related.pairs(x, "parents", available = T)
        p = IBDestimate(x = dat, g1 = P, plot.action = makeplot, pointcol = 2, ...)
    }
    if ("cousins" %in% groups) {
        C = rbind(related.pairs(x, "cousins", half = F, available = T), related.pairs(x, "grandparents", 
            degree = 3, available = T), related.pairs(x, "nephews_nieces", half = T, available = T), 
            related.pairs(x, "nephews_nieces", removal = 2, half = F, available = T))
        co = IBDestimate(x = dat, g1 = C, plot.action = makeplot, pointcol = 6, ...)
    }
    if ("grandparents" %in% groups) {
        G = rbind(related.pairs(x, "siblings", half = T, available = T), related.pairs(x, "grandparents", 
            available = T), related.pairs(x, "nephews_nieces", half = F, available = T))
        g = IBDestimate(x = dat, g1 = G, plot.action = makeplot, pointcol = 3, ...)
    }
    if ("siblings" %in% groups) {
        S = related.pairs(x, "siblings", half = F, available = T)
        s = IBDestimate(x = dat, g1 = S, plot.action = makeplot, pointcol = 4, ...)
    }
    if ("unrelated" %in% groups) {
        interfam = match.arg(interfam)
        U = related.pairs(x, "unrelated", available = T, interfam = interfam)
        u = IBDestimate(x = dat, g1 = U, plot.action = makeplot, pointcol = "darkgray", ...)
    }
    if ("distant" %in% groups) {
        rest = t(combn(x$available, 2))
        taken = rbind(P, U, S, G, C)
        D = rest[!(rest[, 1] * 1000 + rest[, 2]) %in% (taken[, 1] * 1000 + taken[, 2]), , drop = F]
        d = IBDestimate(x = dat, g1 = D, plot.action = makeplot, pointcol = "cyan", ...)
    }
    # if(makeplot) par(oldpar)
    res = list(parents = p, siblings = s, grandparents_uncles_halfsibs = g, cousins = co, distant = d, 
        unrelated = u)
    res = res[!sapply(res, is.null)]
    invisible(res)
}

.IBDsuspects = function(x, radius) {
    # x output from examineKinships, radius=tolerated distance
    
    extract_distant = function(df, x0, y0, rad) {
        if (is.null(df)) 
            return(NULL)
        subs = subset(df, sqrt((df$k0 - x0)^2 + (df$k2 - y0)^2) > rad)
        if (nrow(subs) > 0) 
            subs else NULL
    }
    
    susp = list(`suspects-parents` = extract_distant(x$parents, 0, 0, radius), `suspects-siblings` = extract_distant(x$siblings, 
        0.25, 0.25, radius), `suspects-halfsibs/grandparents/uncles` = extract_distant(x$grandparents, 
        0.5, 0, radius), `suspects-cousins/halfuncles/granduncles` = extract_distant(x$cousins, 
        0.75, 0, radius), `suspects-distant` = extract_distant(x$distant, 0.9, 0, radius), 
        `suspects-unrelated` = extract_distant(x$unrelated, 1, 0, radius))
    susp = susp[!sapply(susp, is.null)]
    if (length(susp) == 0) 
        susp = NULL
}

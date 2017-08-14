#' IBD triangle plot
#' 
#' This function draws an IBD triangle, as typically used to visualize pairwise
#' relatedness of non-inbred individuals.  Various annotations are available,
#' including points marking the most common relationships, contour lines for
#' the kinship coefficients, and shading of the unattainable region.
#' 
#' For any pair of non-inbred individuals A and B, their genetic relationship
#' can be summarized by the IBD coefficients \eqn{(\kappa_0, \kappa_1,
#' \kappa_2)}, where 
#' \deqn{\kappa_i = P(A and B share i alleles IBD at random autosomal locus).}
#' Since \eqn{\kappa_0+\kappa_1+\kappa_2=1}, any
#' relationship corresponds to a point in the triangle in the \eqn{(\kappa_0,
#' \kappa_2)}-plane defined by \eqn{\kappa_0 \ge 0, \kappa_2 \ge 0, \kappa_0 +
#' \kappa_2 \le 1}. The choice of \eqn{\kappa_0} and \eqn{\kappa_2} as the axis
#' variables is done for reasons of symmetry and is not significant (other
#' authors have used different views of the triangle).
#' 
#' As shown in (Thompson, 1976) points in the subset of the triangle defined by
#' \eqn{4\kappa_0\kappa_2 > \kappa_1^2} is unattainable for pairwise
#' relationships.  By default this region in shaded in a 'lightgray' color.
#' 
#' The IBD coefficients are linearly related to the kinship coefficient
#' \eqn{\phi} by the formula \deqn{\phi = 0.25\kappa_1 + 0.5\kappa_2.} By
#' indicating values for \eqn{phi} in the \code{kinship.lines} argument, the
#' corresponding contour lines are shown as dashed lines in the triangle plot.
#' 
#' @param relationships A character vector indicating relationships points to
#' be included in the plot.  UN=unrelated; PO=parent/offspring; MZ=monozygotic
#' twins; S=full siblings; H=half sibling; U=uncle/aunt;
#' @param kinship.lines A numeric vector.
#' @param shading The shading color for the unattainable region.
#' @param pch Symbol used for the relationship points (see \code{\link{par}}).
#' @param cex_points A single numeric controlling the symbol size for the
#' relationship points.
#' @param cex_text A single numeric controlling the font size for the
#' relationship labels.
#' @param axes Draw surrounding axis box?
#' @return NULL
#' @author Magnus Dehli Vigeland
#' @seealso \code{\link{examineKinships}}
#' @references E. A. Thompson (1975). \emph{The estimation of pairwise relationships.}
#' Annals of Human Genetics 39.
#' 
#' E. A. Thompson (1976). \emph{A restriction on the space of genetic relationships.}
#' Annals of Human Genetics 40.
#'
#' @examples
#' 
#' IBDtriangle()
#' 
#' @export
IBDtriangle = function(relationships = c("UN", "PO", "MZ", "S", "H,U,G", "FC", "SC", "DFC", 
    "Q"), kinship.lines = numeric(), shading = "lightgray", pch = 16, cex_points = 1.2, cex_text = 1, 
    axes = FALSE) {
    
    par(xpd = T, mar = c(3.1, 3.1, 1, 1), pty = "s")
    
    plot(NULL, xlim = c(0, 1), ylim = c(0, 1), axes = axes, ann = FALSE)
    mtext(mtext(text = c(expression(italic(kappa[0])), expression(italic(kappa[2]))), side = 1:2, 
        line = c(1, 0.5), las = 1))
    
    # impossible region shading(do borders afterwards)
    kk0 = seq(0, 1, length = 501)
    kk2 = 1 + kk0 - 2 * sqrt(kk0)
    polygon(kk0, kk2, col = shading, border = NA)
    # text(.4, .4, 'impossible region', srt=-45)
    
    # impossible border
    points(kk0, kk2, type = "l", lty = 3)
    
    # axes
    segments(c(0, 0, 0), c(0, 0, 1), c(1, 0, 1), c(0, 1, 0))
    
    # kinship lines
    for (phi in kinship.lines) {
        if (phi < 0 || phi > 0.5) 
            stop(paste("kinship coefficient not in intervall [0, 0.5]", phi))
        abline(a = (4 * phi - 1), b = 1, lty = 2)
        labpos.x = 0.5 * (1.2 - (4 * phi - 1))
        labpos.y = 1.2 - labpos.x
        lab = substitute(paste(phi1, " = ", a), list(a = phi))
        text(labpos.x, labpos.y, labels = lab, pos = 3, srt = 45)
    }
    
    # relationships
    RELS = data.frame(label = c("UN", "PO", "MZ", "S", "H,U,G", "FC", "SC", "DFC", "Q"), k0 = c(1, 
        0, 0, 1/4, 1/2, 3/4, 15/16, 9/16, 17/32), k1 = c(0, 1, 0, 1/2, 1/2, 1/4, 1/16, 6/16, 
        14/32), k2 = c(0, 0, 1, 1/4, 0, 0, 0, 1/16, 1/32), pos = c(1, 1, 4, 4, 1, 1, 1, 3, 2))
    assert_that(is.character(relationships), all(relationships %in% RELS$label))
    if (length(relationships) > 0) {
        rels = RELS[RELS$label %in% relationships, ]
        points(rels$k0, rels$k2, pch = pch, cex = cex_points)
        text(rels$k0, rels$k2, labels = rels$label, pos = rels$pos, cex = cex_text)
    }
}
#' Genotype probability distribution
#'
#' Computes the joint genotype distribution of two markers for a specified
#' pedigree member, conditional on existing genotypes and pedigree information.
#'
#' @param x A \code{\link{linkdat}} object.
#' @param id The individual in question.
#' @param partialmarker1,partialmarker2 Either a single integer indicating the
#'   number of one of \code{x}'s existing markers, or a \code{marker} object.
#' @param theta A single numeric in the interval [0, 0.5] - the recombination
#'   fraction between the two markers.
#' @param loop_breakers A numeric containing IDs of individuals to be used as
#'   loop breakers. Relevant only if the pedigree has loops. See
#'   \code{\link{breakLoops}}.
#' @param eliminate A non-negative integer, indicating the number of iterations
#'   in the internal genotype-compatibility algorithm. Positive values can save
#'   time if \code{partialmarker} is non-empty and the number of alleles is
#'   large.
#' @param verbose A logical.
#' @return A named matrix giving the joint genotype distribution.
#' @seealso \code{\link{oneMarkerDistribution}}
#'
#' @examples
#'
#' x = nuclearPed(2)
#' emptySNP = marker(x, alleles=c('a','b'))
#' SNP1 = marker(x, 1, c(1,1), 2, c(1,0), alleles=1:2, afreq=c(0.1, 0.9))
#' twoMarkerDistribution(x, id=2, emptySNP, SNP1, theta=0)
#' twoMarkerDistribution(x, id=2, emptySNP, SNP1, theta=0.5)
#' twoMarkerDistribution(x, id=3, emptySNP, SNP1, theta=0.5)
#'
#' # X-linked example
#' SNPX_1 = marker(x, 2, c('a','b'), 3, 'b', alleles=c('a','b'), chrom=23)
#' SNPX_2 = marker(x, 2, c('a','b'), 3, 'b', alleles=c('a','b'), chrom=23)
#' r1 = twoMarkerDistribution(x, id=4, SNPX_1, SNPX_2, theta=0)
#' r2 = twoMarkerDistribution(x, id=4, SNPX_1, SNPX_2, theta=0.5)
#' stopifnot(all(r1==c(.5,0,0,.5)), all(r2==c(.25,.25,.25,.25)))
#'
#' @export
twoMarkerDistribution <- function(x, id, partialmarker1, partialmarker2, theta, loop_breakers = NULL,
    eliminate = 99, verbose = TRUE) {
    starttime = proc.time()
    if (!inherits(m1 <- partialmarker1, "marker"))
        if (is.numeric(m1) && length(m1) == 1 && m1 <= x$nMark)
            m1 = x$markerdata[[m1]] else stop("The 'partialmarker1' must be a 'marker' object, or a single integer indicating an existing marker of 'x'.")
    if (!inherits(m2 <- partialmarker2, "marker"))
        if (is.numeric(m2) && length(m2) == 1 && m2 <= x$nMark)
            m2 = x$markerdata[[m2]] else stop("The 'partialmarker2' must be a 'marker' object, or a single integer indicating an existing marker of 'x'.")

    markerchrom1 = as.integer(attr(m1, "chrom"))
    markerchrom2 = as.integer(attr(m2, "chrom"))
    if (!identical(markerchrom1, markerchrom2))
        stop(sprintf("Partial markers are on different chromosomes: %d and %d.", markerchrom1,
            markerchrom2))
    chrom = ifelse(identical(23L, markerchrom1), "X", "AUTOSOMAL")

    SEX = x$pedigree[, "SEX"]
    afreq1 = attr(m1, "afreq")
    alleles1 = attr(m1, "alleles")
    afreq2 = attr(m2, "afreq")
    alleles2 = attr(m2, "alleles")

    if (verbose) {
        cat(ifelse(chrom == "AUTOSOMAL", "Autosomal", "X-linked"), "markers with the following known genotypes:\n")
        print(data.frame(ID = x$orig.ids, M1 = as.character(.prettyMarkers(list(m1), missing = "-",
            singleCol = TRUE, sex = SEX)), M2 = as.character(.prettyMarkers(list(m2), missing = "-",
            singleCol = TRUE, sex = SEX))), row.names = FALSE)
        cat("\nAllele frequencies, marker 1:\n")
        print(structure(afreq1, names = alleles1))
        cat("\nAllele frequencies, marker 2:\n")
        print(structure(afreq2, names = alleles2))
        cat("\nRecombination rate between marker loci:", theta, "\n")
    }

    # Do this before loop breaking, since eliminate2 works better WITH the loops.
    grid.subset = fast.grid(c(geno.grid.subset(x, m1, id, chrom, make.grid = F), geno.grid.subset(x,
        m2, id, chrom, make.grid = F)))

    if (x$hasLoops) {
        x = breakLoops(setMarkers(x, list(m1, m2)), loop_breakers = loop_breakers, verbose = verbose)
        m1 = x$markerdata[[1]]
        m2 = x$markerdata[[2]]
        SEX = x$pedigree[, "SEX"]
    }

    int.id = .internalID(x, id)
    allgenos1 = allGenotypes(attr(m1, "nalleles"))
    allgenos2 = allGenotypes(attr(m2, "nalleles"))

    if (chrom == "AUTOSOMAL" || SEX[int.id] == 2) {
        gt1.strings = paste(alleles1[allgenos1[, 1]], alleles1[allgenos1[, 2]], sep = "/")
        gt2.strings = paste(alleles2[allgenos2[, 1]], alleles2[allgenos2[, 2]], sep = "/")
        geno.names = list(gt1.strings, gt2.strings)
    } else geno.names = list(alleles1, alleles2)

    marginal = likelihood(x, locus1 = m1, locus2 = m2, theta = theta, eliminate = eliminate)
    if (marginal == 0)
        stop("Partial marker data is impossible")

    probs = array(0, dim = lengths(geno.names, use.names = F), dimnames = geno.names)
    probs[grid.subset] = apply(grid.subset, 1, function(allg_rows) {
        m1[int.id, ] = allgenos1[allg_rows[1], ]
        m2[int.id, ] = allgenos2[allg_rows[2], ]
        likelihood(x, locus1 = m1, locus2 = m2, theta = theta, eliminate = eliminate)
    })

    res = probs/marginal
    if (verbose) {
        cat("\nJoint genotype distribution at the two markers for individual ", id, ":\n",
            sep = "")
        print(round(res, 4))
        cat("\nTotal time used: ", (proc.time() - starttime)[["elapsed"]], " seconds.\n", sep = "")
        return(invisible(res))
    } else res
}

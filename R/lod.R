#' Two-point LOD score
#' 
#' Calculates the two-point LOD scores of a pedigree for the specified markers.
#' The recombination ratio between the disease and marker loci can be either
#' fixed at specific values, or optimized.
#' 
#' The LOD score of a marker is defined as 
#' \deqn{LOD(\theta) = \log[10] \frac{L(\theta)}{L(0.5)}}{%
#' LOD(\theta) = log10 L(\theta)/L(0.5)}
#' where \eqn{L(\theta)} denotes the likelihood of the
#' observed marker genotypes given a recombination ratio \eqn{\theta}
#' between the marker and the disease locus.
#' 
#' @param x a \code{\link{linkdat}} object.
#' @param markers an integer vector denoting which markers to use.
#' @param theta either a numeric containing specific recombination ratio(s), or
#' the word 'max', indicating that the recombination ratio should be optimized
#' by the program.
#' @param loop_breakers a numeric containing IDs of individuals to be used as
#' loop breakers. Relevant only if the pedigree has loops. See
#' \code{\link{breakLoops}}.
#' @param max.only a logical indicating whether only the maximum LOD score
#' should be returned.
#' @param verbose a logical: verbose output or not.
#' @param tol a numeric passed on to \code{\link{optimize}} as its tolerance
#' parameter.
#'
#' @return If \code{max.only=TRUE}, the highest computed LOD score is returned,
#'  as a single number.
#'
#' Otherwise a \code{linkres} object, which is essentially a matrix containing 
#' the LOD scores. The details depend on the other parameters:
#'
#' If \code{theta} is numeric, the matrix has
#' dimensions \code{length(theta) * length(markers)}, and the entry in row
#' \code{t}, column \code{m} is the lod score of the pedigree for marker
#' \code{m} assuming a recombination rate of \code{t}.
#' 
#' If \code{theta='max'}, the \code{linkres} matrix has one column per marker
#' and two rows: The first containing the LOD score and the second the optimal
#' recombination ratio for each marker.
#' 
#' If a marker has incompatible values (i.e. if a child of homozygous 1/1
#' parents has a 2 allele), the corresponding output entries are \code{NaN}.
#' 
#' @author Magnus Dehli Vigeland
#' @seealso \code{\link{likelihood}}, \code{\link{optimize}},
#' \code{\link{breakLoops}}
#' 
#' @examples
#' 
#' x = linkdat(toyped, model=1)
#' res = lod(x)
#' res_theta = lod(x, theta=c(0, 0.1, 0.2, 0.5))
#' res_max = lod(x, theta='max')
#' stopifnot(all(0.3 == round(c(res, res_theta['0',], res_max['LOD',]), 1)))
#' 
#' # bigger pedigree with several markers
#' y = linkdat(dominant)
#' y = setModel(y, model=1, penetrances=c(.001, .9, .99))
#' lod(y, markers=305:310)
#' lod(y, markers=305:310, theta='max')
#' 
#' # Example with pedigree with loops:
#' z = linkdat(twoloops, model=2) # fully penetrant autosomal recessive model.
#' 
#' # add SNP for which individuals 15 and 16 are homozygous for the rare allele.
#' m = marker(z, 15:16, c(1,1), alleles=1:2, afreq=c(0.001, 0.999))
#' z = addMarker(z, m)
#' res1 = lod(z)
#' # manual specification of loop breakers gives same result
#' res2 = lod(z, loop_breakers=c(8,12))
#' 
#' # making the marker triallelic and adding some genotypes.
#' z = modifyMarker(z, marker=1, ids=c(7,9,11,13), genotype=3, alleles=1:3, afreq=c(0.001, 0.499, 0.5))
#' plot(z, 1)
#' res3 = lod(z)
#' 
#' z = modifyMarker(z, marker=1, alleles=1:4, afreq=c(0.001, 0.499, 0.25, 0.25)) 
#' res4 = lod(z)
#' 
#' stopifnot(all(3 == round(c(res1, res2, res3, res4), 1)))
#' 
#' @export
lod = function(x, markers = seq_len(x$nMark), theta = 0, loop_breakers = NULL, max.only = FALSE, 
    verbose = FALSE, tol = 0.01) {
    assert_that(is.linkdat(x))
    if (is.singleton(x)) 
        stop("This function is not applicable to singleton objects.")
    if (is.null(x$model)) 
        stop("No model set.")
    if (x$nMark == 0) 
        stop("No marker data indicated.")
    if (x$hasLoops) {
        x = breakLoops(x, loop_breakers = loop_breakers, verbose = verbose)
    }
    if (min(markers) < 1 || max(markers) > x$nMark) 
        stop("Nonexistent marker indicated.")
    chrom = x$model$chrom
    
    maxall = max(unlist(lapply(x$markerdata[markers], attr, "nalleles")))
    if (maxall > 4) {
        if (is.character(theta) && theta == "max") 
            stop("The option theta='max' is not implemented for markers with more than 4 alleles.")
        return(.lodm(x, markers = markers, theta = theta, loop_breakers = loop_breakers, max.only = max.only, 
            verbose = verbose))
    }
    
    ilink <- function(x, marker) {
        m = x$markerdata[[marker]]
        nall = attr(m, "nalleles")
        afreq = attr(m, "afreq")
        initialCalc = .initialCalc(x, afreq = afreq, chrom)
        
        log_denom = stopl = likelihood_LINKAGE(x, m, logbase = 10, initialCalc = initialCalc, 
            TR.MATR = TR05[[nall]])
        if (log_denom == -Inf) 
            return(c(NaN, NaN))
        
        startl = likelihood_LINKAGE(x, m, logbase = 10, TR.MATR = TRZ[[nall]])
        optimal = optimize(likelihood_LINKAGE, c(0, 0.5), x = x, marker = m, afreq = NULL, 
            initialCalc = initialCalc, logbase = 10, tol = tol, maximum = TRUE)
        if (optimal$objective > max(startl, stopl)) {
            log_numer <- optimal$objective
            theta_max <- optimal$maximum
        } else {
            log_numer = max(startl, stopl)
            theta_max <- c(0, 0.5)[which.max(c(startl, stopl))]
        }
        c(log_numer - log_denom, theta_max)
    }
    
    TR05 = lapply(1:maxall, function(n) .TRmatrNEW(0.5, n, chrom))
    
    map = .getMap(x, na.action = 1, verbose = F)[markers, , drop = F]
    
    if (is.numeric(theta)) {
        stopifnot(max(theta) <= 0.5, min(theta) >= 0)
        if (verbose) 
            cat(sprintf("Computing singlepoint LOD scores at each marker\nfor the following recombination value%s:\n%s\n", 
                ifelse(length(theta) == 1, "", "s"), paste("  theta =", theta, collapse = "\n")))
        
        markerdata_list = x$markerdata[markers]
        trm_list = lapply(theta, function(th) lapply(1:maxall, function(n) .TRmatrNEW(th, n, 
            chrom)))
        denoms = unlist(lapply(markerdata_list, function(m) likelihood_LINKAGE(x, m, theta = NULL, 
            logbase = 10, TR.MATR = TR05[[attr(m, "nalleles")]])))
        numers = vapply(markerdata_list, function(m) unlist(lapply(trm_list, function(TR) likelihood_LINKAGE(x, 
            marker = m, theta = NULL, logbase = 10, TR.MATR = TR[[attr(m, "nalleles")]]))), 
            FUN.VALUE = numeric(length(theta)))
        res = numers - rep(denoms, each = length(theta))
        res = structure(res, dim = c(length(theta), length(markers)), dimnames = list(theta, 
            map$MARKER), analysis = "mlink", map = map, class = "linkres")
    } else if (identical(theta, "max")) {
        if (verbose) 
            cat("Computing singlepoint LOD scores for each marker,\nmaximizing over all recombination values.\n\n")
        TRZ = lapply(1:maxall, function(n) .TRmatrNEW(0, n, chrom))
        res = unlist(lapply(markers, ilink, x = x))
        res = structure(res, dim = c(2, length(markers)), dimnames = list(c("LOD", "t_max"), 
            map$MARKER), analysis = "ilink", map = map, class = "linkres")
    }
    
    if (verbose) 
        summary.linkres(res)
    if (max.only) 
        ifelse(all(is.na(res)), NA, max(res, na.rm = TRUE)) else res
}




.lodm = function(x, markers = seq_len(x$nMark), theta = 0, loop_breakers = NULL, max.only = FALSE, 
    verbose = FALSE) {
    assert_that(is.linkdat(x))
    if (is.null(x$model)) 
        stop("No model set.")
    if (x$nMark == 0) 
        stop("No marker data indicated.")
    if (x$hasLoops) {
        x = breakLoops(x, loop_breakers = loop_breakers, verbose = verbose)
    }
    
    theta_list = as.list(theta)
    lods = vapply(x$markerdata[markers], function(marker) {
        attr(marker, "chrom") = x$model$chrom
        start.dat = .startdata_MD(x, marker)
        denom = likelihood.linkdat(x, locus1 = marker, locus2 = "disease", theta = 0.5, startdata = start.dat, 
            logbase = 10)
        unlist(lapply(theta_list, function(tht) likelihood.linkdat(x, locus1 = marker, locus2 = "disease", 
            theta = tht, startdata = start.dat, logbase = 10) - denom))
    }, FUN.VALUE = numeric(length(theta)))
    
    map = .getMap(x, na.action = 1, verbose = F)[markers, , drop = FALSE]
    
    res = structure(lods, dim = c(length(theta), length(markers)), dimnames = list(theta, map$MARKER), 
        analysis = "mlink", map = map, class = "linkres")
    if (verbose) 
        summary(res)
    if (max.only) 
        ifelse(all(is.na(res)), NA, max(res, na.rm = TRUE)) else res
}


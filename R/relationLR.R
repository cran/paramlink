#' Relationship Likelihood Ratio
#' 
#' Computes likelihood for two pedigrees and their ratio, the likelihood ratio
#' (LR).
#' 
#' This function computes the likelihood of two pedigrees (each corresponding
#' to a hypothesis describing a family relationship).  The likelihood ratio is
#' also reported. Unlike other implementations we are aware of, partial DNA
#' profiles are allowed here. For instance, if the genotype of a person is
#' reported as 1/0 (0 is 'missing') for a triallelic marker with uniform allele
#' frequencies, the possible ordered genotypes (1,1), (2,1), (1,2), (1,3) and
#' (3,1) are treated as equally likely. (For general allele frequencies,
#' genotype probabilities are obtained by assuming Hardy-Weinberg equilibrium.)
#' A reasonable future extension would be to allow the user to weigh these
#' genotypes; typically (1,1) may be more likely than the others.  If
#' \code{plot='plot_only'}, the function returns NULL after producing the plot.
#' 
#' @param ped_numerator a \code{\link{linkdat}} object, or a list of several
#' linkdat and/or singleton objects, describing the relationship corresponding
#' to the hypothesis H1 (numerator).  If a list, the sets of ID labels must be
#' disjoint, that is, all ID labels must be unique.
#' @param ped_denominator a \code{\link{linkdat}} object, or a list of several
#' linkdat and/or singleton objects, describing the relationship corresponding
#' to the hypothesis H2 (denominator). ID labels must be consistent with
#' \code{ped_claim}.
#' @param ids genotyped individuals.
#' @param alleles a numeric or character vector containing marker alleles names
#' @param afreq a numerical vector with allele frequencies. An error is given
#' if they don't sum to 1 (rounded to 3 decimals).
#' @param known_genotypes list of triplets \code{(a, b, c)}, indicating that
#' individual \code{a} has genotype \code{b/c}. Missing value is 0.
#' @param loop_breakers Not yet implemented, only default value NULL currently
#' handled
#' @param Xchrom a logical: Is the marker on the X chromosome?
#' @param plot either a logical or the character 'plot_only', controlling if a
#' plot should be produced. If 'plot_only', a plot is drawn, but no further
#' computations are done.
#' @param title1 a character, title of leftmost plot.
#' @param title2 a character, title of rightmost plot.
#' @return \item{lik.numerator}{likelihood of data given ped_numerator}
#' \item{lik.denominator}{likelihood of data given ped_denominator}
#' \item{LR}{likelihood ratio lik.numerator/lik.denominator}
#' @author Thore Egeland, Magnus Dehli Vigeland
#' @seealso \code{\link{exclusionPower}}
#' 
#' @examples
#' 
#' ############################################
#' # A partial DNA profile is obtained from the person 
#' # denoted 4 in the figure produced below 
#' # There are two possibilities: 
#' # H1: 4 is the missing relative of 3 and 6 (as shown to the left) or 
#' # H2: 4 is unrelated to 3 and 6.
#' ############################################
#' p = c(0.2, 0.8)
#' alleles = 1:length(p)
#' g3 = c(1,1); g4 = c(1,0); g6 = c(2,2)
#' x1 = nuclearPed(2)
#' x1 = addOffspring(x1, father = 4, sex = 1, noff = 1)
#' m = marker(x1, 3, g3, 4, g4, 6, g6, alleles = alleles, afreq = p)
#' x1 = addMarker(x1, m)
#' x2 = nuclearPed(2)
#' x2 = addOffspring(x2, father = 4, sex = 1, noff = 1)
#' m = marker(x2, 3, g3, 6, g6, alleles = alleles, afreq = p)
#' x2 = addMarker(x2, m)
#' missing = singleton(4, sex = 1)
#' m.miss = marker(missing, g4, alleles = alleles, afreq = p)
#' missing = addMarker(missing, m.miss)
#' x2 = relabel(x2, c(1:3, 99, 5:6), 1:6)
#' known = list(c(3, g3), c(4,g4), c(6, g6))
#' LR = relationLR(x1, list(x2, missing), ids = c(3,4,6),
#'                 alleles = alleles, afreq = p, known = known,
#'                 title1 = 'H1: Missing person 4 related',
#'                 title2 = 'H2:Missing person 4 unrelated')$LR
#' # Formula:
#' p = p[1]
#' LR.a = (1+p)/(2*p*(2-p))
#' stopifnot(abs(LR - LR.a) < 1e-10)
#' 
#' @export
relationLR = function(ped_numerator, ped_denominator, ids, alleles, afreq = NULL, known_genotypes = list(), 
    loop_breakers = NULL, Xchrom = FALSE, plot = TRUE, title1 = "", title2 = "") {
    
    ped_claim = ped_numerator
    ped_true = ped_denominator
    if (is.linkdat(ped_claim)) 
        ped_claim = list(ped_claim)
    if (is.linkdat(ped_true)) 
        ped_true = list(ped_true)
    ids_claim = lapply(ped_claim, function(x) ids[ids %in% x$orig.ids])
    ids_true = lapply(ped_true, function(x) ids[ids %in% x$orig.ids])
    
    if (!is.null(loop_breakers)) 
        return("Loops not yet implemented")
    # loops_claim = lapply(ped_claim, function(x) { lb = x$orig.ids[x$orig.ids %in%
    # loop_breakers] if (length(lb) == 0) lb = NULL lb }) loops_true = lapply(ped_true,
    # function(x) { lb = x$orig.ids[x$orig.ids %in% loop_breakers] if (length(lb) == 0) lb =
    # NULL lb })
    
    N_claim = length(ped_claim)
    N_true = length(ped_true)
    N = N_claim + N_true
    if (length(alleles) == 1) 
        alleles = seq_len(alleles)
    chrom = if (Xchrom) 
        23 else NA
    partial_claim = lapply(1:N_claim, function(i) {
        x = ped_claim[[i]]
        m = marker(x, alleles = alleles, afreq = afreq, chrom = chrom)
        for (tup in known_genotypes) if (tup[1] %in% x$orig.ids) 
            m = modifyMarker(x, m, ids = tup[1], genotype = tup[2:3])
        m
    })
    partial_true = lapply(1:N_true, function(i) {
        x = ped_true[[i]]
        m = marker(x, alleles = alleles, afreq = afreq, chrom = chrom)
        for (tup in known_genotypes) if (tup[1] %in% x$orig.ids) 
            m = modifyMarker(x, m, ids = tup[1], genotype = tup[2:3])
        m
    })
    if (isTRUE(plot) || plot == "plot_only") {
        op = par(oma = c(0, 0, 3, 0), xpd = NA)
        widths = ifelse(sapply(c(ped_claim, ped_true), is.singleton), 1, 2)
        claim_ratio = sum(widths[1:N_claim])/sum(widths)
        layout(rbind(1:N), widths = widths)
        has_genotypes = length(known_genotypes) > 0
        for (i in 1:N) {
            if (i <= N_claim) {
                x = ped_claim[[i]]
                avail = ids_claim[[i]]
                mm = if (has_genotypes) 
                  partial_claim[[i]] else NULL
            } else {
                x = ped_true[[i - N_claim]]
                avail = ids_true[[i - N_claim]]
                mm = if (has_genotypes) 
                  partial_true[[i - N_claim]] else NULL
            }
            cols = ifelse(x$orig.ids %in% avail, 2, 1)
            plot(x, marker = mm, col = cols, margin = c(2, 4, 2, 4), title = "")
        }
        mtext(title1, outer = TRUE, at = claim_ratio/2)
        mtext(title2, outer = TRUE, at = 0.5 + claim_ratio/2)
        rect(grconvertX(0.02, from = "ndc"), grconvertY(0.02, from = "ndc"), grconvertX(claim_ratio - 
            0.02, from = "ndc"), grconvertY(0.98, from = "ndc"))
        rect(grconvertX(claim_ratio + 0.02, from = "ndc"), grconvertY(0.02, from = "ndc"), 
            grconvertX(0.98, from = "ndc"), grconvertY(0.98, from = "ndc"))
        par(op)
        if (plot == "plot_only") 
            return()
    }
    
    lik.claim = lik.true = 1
    for (i in 1:N_claim) lik.claim = lik.claim * likelihood(ped_claim[[i]], partial_claim[[i]])
    for (i in 1:N_true) lik.true = lik.true * likelihood(ped_true[[i]], partial_true[[i]])
    list(lik.numerator = lik.claim, lik.denominator = lik.true, LR = lik.claim/lik.true)
}

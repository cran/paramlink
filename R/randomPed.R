#' Random pedigree
#' 
#' Creates a random medical pedigree with specified number of generations.
#' 
#' The function produces a random simple pedigree. Each founder is given at
#' most one disease allele. At least one of the two top founders carries a
#' disease allele.
#' 
#' @param gen an integer in the interval \code{[2,5]} indicating the number of
#' generations.
#' @param lambda a positive numeric. For each descendant of the first
#' generation, the number of offspring is sampled from a Poisson distribution
#' with parameter \code{lambda}.
#' @param penetrances a numeric of length 3, as in \code{\link{setModel}}.
#' @param naff an integer specifying a lower bound on the number of affected
#' individuals, or the character 'last.gen'. The latter produce a pedigree
#' where at least one in the youngest generation is affected.
#' @param founder.mut an integer, the number of disease alleles to be
#' distributed among the founders.
#' @return A \code{linkdat} object.
#' @author Magnus Dehli Vigeland
#' @seealso \code{\link{linkdat}}
#' 
#' @examples
#' 
#' plot(randomPed(3))
#' 
#' # gives error message: Not enough founder mutations
#' \dontrun{
#' randomPed(gen=4, penetrances=c(0,0,1), naff=2, founder.mut=1)
#' }
#' 
#' @export
randomPed = function(gen, lambda = 2, penetrances = c(0, 1, 1), naff = "last.gen", founder.mut = 1) {
    stopifnot(gen %in% 2:5, length(penetrances) == 3, any(penetrances > 0))
    if (sum(penetrances[1:2]) == 0 && founder.mut < 2) 
        stop("Recessive disease: More founder mutations needed")
    add.offsp <- function(ped, id, offs) {
        if (offs == 0) 
            return(ped)
        N = nrow(ped)
        G = ped[[id, "GENERATION"]]
        S = ped[id, "SEX"]
        parents = c(id, N + 1)[order(c(S, 3 - S))]
        rbind(ped, c(N + 1, 0, 0, 3 - S, 1, G), cbind(ID = (N + 2):(N + offs + 1), FID = parents[1], 
            MID = parents[2], SEX = sample.int(2, offs, replace = TRUE), AFF = 1, GENERATION = G + 
                1))
    }
    
    #-----simulate pedigree---------
    rpois.mix = function(n, lamb, p = 0.5) ifelse(runif(n) < p, 0, rpois(n, lamb))  #mixed fordeling: 50% poisson(lambda)
    
    founder = matrix(c(1, 0, 0, 1, 1, 1), nrow = 1, dimnames = list(NULL, c("ID", "FID", "MID", 
        "SEX", "AFF", "GENERATION")))
    ped = add.offsp(founder, 1, rpois(1, lambda) + 1)  #ensure second generation exists
    
    while ((g <- ped[[nrow(ped), "GENERATION"]]) < gen) for (f in which(ped[, "GENERATION"] == 
        g)) ped = add.offsp(ped, f, offs = rpois.mix(n = 1, lambda))
    
    #------distribute mutations and simulate affection status---------
    founders <- which(ped[, "FID"] == 0)
    nonfounders <- which(ped[, "FID"] > 0)
    cond = FALSE
    tries = 0
    maxtries = 100
    while (!cond && (tries <- (tries + 1)) < maxtries) {
        alleles = matrix(1, nrow = nrow(ped), ncol = 2)
        carrier = sample.int(2, 1)
        if (founder.mut > 1) {
            # If specified, select additional founder mutation carriers
            carrier = c(carrier, sample(setdiff(founders, carrier), size = min(founder.mut, 
                length(founders)) - 1))
        }
        alleles[carrier, 1] <- 2
        for (i in nonfounders) {
            fid = ped[i, "FID"]
            mid = ped[i, "MID"]  # Identify parents
            alleles[i, 1] <- sample(alleles[fid, ], size = 1)  # Draw allele from father
            alleles[i, 2] <- sample(alleles[mid, ], size = 1)  # Draw allele from mother
        }
        dis_alleles = rowSums(alleles - 1)  # Number of disease alleles for each indiv 
        ped[, "AFF"] = sapply(dis_alleles, function(d) sample.int(2, size = 1, prob = c(1 - 
            penetrances[d + 1], penetrances[d + 1])))
        
        if (naff == "last.gen") 
            cond = any(as.logical(ped[ped[, "GENERATION"] == gen, "AFF"] - 1))  #check whether last generation is affected
 else if (is.numeric(naff)) 
            cond = (sum(ped[, "AFF"] - 1) >= naff)  #check whether at least 'naff' are affected
    }
    if (tries == maxtries) 
        return(randomPed(gen = gen, lambda = lambda, penetrances = penetrances, naff = naff, 
            founder.mut = founder.mut)) else return(linkdat(ped[, 1:5], verbose = FALSE))
}

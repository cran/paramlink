#' Likelihood ratios of pedigree hypotheses
#'
#' This function computes likelihood ratios for a given a list of pedigrees (linkdat/singletons objects), one of which is
#' the 'reference', with genotype data from the same set of markers. Data exported from the 'Familias' software can be analysed
#' by using \code{\link{Familias2linkdat}} prior to calling this function.
#' 
#' @param x A list of pedigrees. Each pedigree is either a single linkdat/singleton object, or a list of such objects 
#' (the latter is neccessary if the pedigree is disconnected).
#' @param ref A single integer, indicating the index of the reference pedigree. This is used in the denominator of each LR.
#' @param markers A vector of integers, indexing which markers should be included. If NULL (the default) all markers are used.
#' @return A list with entries
#' \item{LR}{Likelihood ratios}
#' \item{LRperMarker}{Likelihood ratios for each marker}
#' \item{likelihoodsPerSystem}{Likelihoods for each marker}
#' \item{time}{user, system and elapsed time}
#' @author Magnus Dehli Vigeland and Thore Egeland 
#' @seealso \code{\link{IBDtriangle}}, \code{\link{examineKinships}}
#'
#' @examples
#' 
#' # Simulate genotypes for 5 tetraallelic markers for a pair of full sibs
#' set.seed(123)
#' sibs = simpleSim(nuclearPed(2), N=5, alleles=1:4, available=3:4)
#' 
#' # Create two alternative hypotheses and transfer the simulated genotypes to them
#' halfsibs = addOffspring(nuclearPed(1),father=1,noffs=1,id=4)
#' halfsibs = transferMarkerdata(sibs, halfsibs)
#'
#' unrel = list(singleton(3), singleton(4))
#' unrel = transferMarkerdata(sibs, unrel)
#'
#' # Compute LR with 'unrelated' as reference
#' LR(list(sibs, halfsibs, unrel), ref=3)
#'
#' \dontrun{
#' data(adoption)
#' x = Familias2linkdat(adoption$pedigrees, adoption$datamatrix, adoption$loci)
#' result = LRparamlink(x, ref=2)
#' result33 = LRparamlink(x, ref=2, marker=c(11,33))
#' }
#' 
#' @export
LR = function(x, ref, markers) {
    st = proc.time()
  
    # get marker names
    y = if(is.linkdat(x[[1]])) x[[1]] else x[[1]][[1]]
    if(missing(markers)) 
        markers = seq_len(y$nMark)
    markernames = sapply(y$markerdata[markers], attr, "name")
    NAnames = is.na(markernames)
    if(any(NAnames)) markernames[NAnames] = paste0("M", which(NAnames))
    
    # break all loops (NB: would use rapply, but doesnt work since is.list(linkdat) = TRUE
    breaklp = function(a) 
        if(class(a)[[1]] == "linkdat") # not singletons
            breakLoops(a, verbose=F) 
        else 
            a 

    x_loopfree = lapply(x, function(xx) if(is.linkdat.list(xx))  lapply(xx, breaklp) else breaklp(xx))

    # compute likelihoods
    liks = lapply(x_loopfree, function(xx) vapply(markers, function(i)
        likelihood(xx, locus1=i), FUN.VALUE=1))
    likelihoodsPerSystem = do.call(cbind, liks)

    # LR per marker and total
    LRperMarker = do.call(cbind, 
        lapply(1:length(x), function(j) liks[[j]]/liks[[ref]]))
    
    # total LR
    LR = apply(LRperMarker, 2, prod)

    # output
    time = proc.time()-st
    rownames(likelihoodsPerSystem) = rownames(LRperMarker) = markernames
    list(LR=LR, LRperMarker=LRperMarker, likelihoodsPerSystem=likelihoodsPerSystem, time=time)
}


#' Check for Mendelian errors
#'
#' Check marker data for Mendelian inconsistencies
#'
#' @param x a \code{\link{linkdat}} object
#' @param remove a logical. If FALSE, the function returns the indices of
#' markers found to incorrect.  If TRUE, a new \code{linkdat} object is
#' returned, where the incorrect markers have been deleted.
#' @param verbose a logical. If TRUE, details of the markers failing the tests
#' are shown.
#'
#' @return
#' A numeric containing the indices of the markers
#' that did not pass the tests, or (if \code{remove=TRUE}) a new \code{linkdat}
#' object where the failing markers are removed.
#' @author Magnus Dehli Vigeland
#' @examples
#' 
#' x = nuclearPed(3)
#' 
#' # Adding a SNP with a mendelian error:
#' # Individual 3 has an allele 'c' not carried by either parents
#' m1 = marker(x, 1, c('a','a'), 2, c('a','b'), 3, c('a','c')) 
#' 
#' # Another erroneous marker: The siblings carry more than 4 different alleles.
#' m2 = marker(x, 3, c(1,2), 4, c(3,4), 5, c(1,5))
#' 
#' # Another marker with incosistent genotypes among the siblings:
#' m3 = marker(x, 3, c(1,1), 4, c(2,2), 5, c(3,3))
#' 
#' # Another marker with incosistent genotypes among the siblings:
#' m4 = marker(x, 3, c(1,1), 4, c(2,3), 5, c(1,4))
#' 
#' # A correct marker (all homozygous for allele 'A')
#' m5 = marker(x, 1:5, 'A')
#'
#' # An empty marker
#' m6 = marker(x)
#'
#' x = setMarkers(x, list(m1,m2,m3,m4,m5,m6))
#' 
#' # Finding the errors
#' err_index = mendelianCheck(x, remove=FALSE)
#' stopifnot(all.equal(err_index, 1:4))
#' 
#' x_remove = mendelianCheck(x, remove=TRUE)
#' stopifnot(x_remove$nMark == 2)
#' 
#' @export
mendelianCheck = function(x, remove = FALSE, verbose = !remove) {
    if (is.singleton(x)) 
        return(numeric(0))
    
    trioCheckFast = function(fa, mo, of) {
        even = 2 * seq_len(length(fa)/2)
        odd = even - 1
        fa_odd = fa[odd]
        fa_even = fa[even]
        mo_odd = mo[odd]
        mo_even = mo[even]
        of_odd = of[odd]
        of_even = of[even]
        fa0 = (fa_odd == 0 | fa_even == 0)
        mo0 = (mo_odd == 0 | mo_even == 0)
        of_odd0 = (of_odd == 0)
        of_even0 = (of_even == 0)
        ff1 = (fa0 | of_odd0 | of_odd == fa_odd | of_odd == fa_even)
        ff2 = (fa0 | of_even0 | of_even == fa_odd | of_even == fa_even)
        mm1 = (mo0 | of_odd0 | of_odd == mo_odd | of_odd == mo_even)
        mm2 = (mo0 | of_even0 | of_even == mo_odd | of_even == mo_even)
        (ff1 & mm2) | (ff2 & mm1)
    }
    
    maleXHomoz = function(of) {
        even = 2 * seq_len(length(of)/2)
        odd = even - 1
        of[odd] == of[even]
    }
    
    maleXCheck = function(mo, of) {
        even = 2 * seq_len(length(of)/2)
        odd = even - 1
        mo_odd = mo[odd]
        mo_even = mo[even]
        of = of[odd]
        mo0 = (mo_odd == 0 | mo_even == 0)
        of == 0 | mo0 | of == mo_odd | of == mo_even
    }
    
    sibshipCheck = function(offs) {
        # offs = matrix with 2*N columns
        even = 2 * seq_len(ncol(offs)/2)
        
        # loop through markers
        unlist(lapply(even, function(i) {
            ###offs_als = unique.default(offs[, (i - 1):i])
            genos = offs[, (i - 1):i]
            
            # number of (different) alleles occuring in homozygous state
            homoz_alleles = genos[genos[,1] == genos[,2], 1]
            n_homoz = length(.mysetdiff(homoz_alleles, 0))
            
            # number of different alleles in total
            n_alleles = length(.mysetdiff(genos, 0))
            
            # if no homoz: consistent if number of alleles <= 4.
            # for each new allele observed homoz, "4" is reduced by 1.
            n_alleles <= (4 - n_homoz)
        }))
    }
    
    ped = x$pedigree
    parents = unique(ped[, 2:3])
    parents = parents[-match(0, parents[, 1]), , drop = FALSE]
    subnucs = lapply(nrow(parents):1, function(i) {
        par = parents[i, ]
        c(fa = par[[1]], mo = par[[2]], offs = as.vector(ped[, 1])[which(ped[, 2] == par[[1]] & 
            ped[, 3] == par[[2]], useNames = FALSE)])
    })
    
    chromX = unlist(lapply(x$markerdata, function(mm) identical(23L, as.integer(attr(mm, "chrom")))))
    which_AUT = which(!chromX)
    which_X = which(chromX)
    
    errorlist = rep(list(numeric(0)), nrow(ped))
    names(errorlist) = x$orig.ids
    nuc_errors = numeric()  # container for allele count errors...belongs to the whole subnuc.
    
    ### AUTOSOMAL
    if (length(which_AUT) > 0) {
        if (verbose) 
            cat("\n### Checking autosomal markers ###\n")
        mdat = do.call(cbind, x$markerdata[!chromX])
        for (sub in subnucs) {
            fa = mdat[sub[[1]], ]
            mo = mdat[sub[[2]], ]
            offs = sub[-(1:2)]
            for (of in offs) {
                new_errors = which_AUT[!trioCheckFast(fa, mo, mdat[of, ])]
                if (length(new_errors) > 0) {
                  errorlist[[of]] = c(errorlist[[of]], new_errors)
                  if (verbose) 
                    cat(sprintf("Individual %d incompatible with parents for %d markers: %s\n", 
                      x$orig.ids[of], length(new_errors), paste(new_errors, collapse = ", ")))
                }
            }
            if (length(offs) > 1) {
                new_errors = which_AUT[!sibshipCheck(offs = mdat[sub[-(1:2)], ])]
                if (length(new_errors) > 0) {
                  nuc_errors = c(nuc_errors, new_errors)
                  if (verbose) 
                    cat(sprintf("Offspring of %d and %d have too many alleles for %d markers: %s\n", 
                      x$orig.ids[sub[1]], x$orig.ids[sub[2]], length(new_errors), paste(new_errors, 
                        collapse = ", ")))
                }
            }
        }
    }
    
    ### X
    if (length(which_X) > 0) {
        if (verbose) 
            cat("\n### Checking markers on the X chromosome ###\n")
        sex = x$pedigree[, "SEX"]
        mdat = do.call(cbind, x$markerdata[chromX])
        for (sub in subnucs) {
            fa = mdat[sub[[1]], ]
            mo = mdat[sub[[2]], ]
            offs = sub[-(1:2)]
            for (of in offs) {
                ofdat = mdat[of, ]
                if (sex[of] == 1) {
                  new_errors = which_X[!maleXHomoz(ofdat)]
                  if (length(new_errors) > 0) {
                    errorlist[[of]] = c(errorlist[[of]], new_errors)
                    if (verbose) 
                      cat("Male ", x$orig.ids[of], " heterozygous for ", length(new_errors), 
                        " X-linked markers: ", paste(new_errors, collapse = ", "), "\n", sep = "")
                  }
                  new_errors = which_X[!maleXCheck(mo, ofdat)]
                  if (length(new_errors) > 0) {
                    errorlist[[of]] = c(errorlist[[of]], new_errors)
                    if (verbose) 
                      cat("Male ", x$orig.ids[of], " incompatible with mother for ", length(new_errors), 
                        " X-linked markers: ", paste(new_errors, collapse = ", "), "\n", sep = "")
                  }
                } else {
                  new_errors = which_X[!trioCheckFast(fa, mo, ofdat)]
                  if (length(new_errors) > 0) {
                    errorlist[[of]] = c(errorlist[[of]], new_errors)
                    if (verbose) 
                      cat("Female ", x$orig.ids[of], " incompatible with parents for ", length(new_errors), 
                        " X-linked markers: ", paste(new_errors, collapse = ", "), "\n", sep = "")
                  }
                }
            }
            if (length(offs) > 2) {
                new_errors = which_X[!sibshipCheck(offs = mdat[offs, ])]  #TODO: fix this for X-linked
                if (length(new_errors) > 0) {
                  nuc_errors = c(nuc_errors, new_errors)
                  if (verbose) 
                    cat(sprintf("Offspring of %d and %d have too many alleles for %d markers: %s\n", 
                      x$orig.ids[sub[1]], x$orig.ids[sub[2]], length(new_errors), paste(new_errors, 
                        collapse = ", ")))
                  
                }
            }
        }
    }
    err_index = sort.int(unique.default(c(unlist(errorlist), nuc_errors)))
    
    if (remove) 
        return(removeMarkers(x, err_index)) else return(err_index)
}

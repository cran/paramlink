#' Relatedness estimation
#'
#' Estimate the pairwise IBD coefficients \eqn{(\kappa_0, \kappa_1, \kappa_2)} for
#' specified pairs of pedigree members, using maximum likelihood methods.
#' The optimization machinery is imported from the \code{maxLik} package.
#'
#' This function optimises the log-likelihood function first described in (Thompson, 1975).
#' Optimisation is done in the \eqn{(\kappa_0, \kappa_2)}-plane and restricted to the
#' probability triangle defined by \eqn{\kappa_0 \ge 0, \kappa_2 \ge 0, \kappa_0 + \kappa_2 \le 1}.

#' @param x A single \code{linkdat} object or a list of \code{linkdat} and/or
#'  \code{singleton} objects.
#' @param ids Either a vector of length 2, or a matrix with two columns, indicating the
#'  the pair(s) of individuals for which IBD estimates should be computed. If a matrix,
#'  each row corresponds to a pair. The entries can be either characters (matching the
#'  \code{plot.labels} of the linkdat object(s)) or integers (matching the \code{orig.ids}
#'  identifiers of the linkdat object(s)).
#' @param markers A numeric indicating which marker(s) to include. If NULL (default),
#'  all markers are used.
#' @param start Numeric of length 2, indicating the initial value of \eqn{(\kappa_0, \kappa_2)}
#'  in the optimisation (passed on to \code{maxLik}).
#' @param tol A single numeric: the optimising tolerance value; passed on to \code{maxLik}).
#'
#' @return A data.frame with 8 columns: ID1, ID2 (numeric IDs), Name1, Name2 (plot labels, if present), N
#' (#markers with no missing alleles), \eqn{\kappa_0}, \eqn{\kappa_1}, \eqn{\kappa_2}.
#' @seealso \code{\link{examineKinships}},
#'  \code{\link{IBDtriangle}}, \code{\link[maxLik]{maxLik}}
#'
#' @references E. A. Thompson (2000). \emph{Statistical Inferences from Genetic Data
#'  on Pedigrees.} NSF-CBMS Regional Conference Series in Probability and
#'  Statistics. Volume 6.
#'
#' @examples
#'
#' if (requireNamespace("maxLik", quietly = TRUE)) {
#'
#' # Simulate marker data for two siblings
#' x = nuclearPed(2)
#' x = setPlotLabels(x, c("Sib1", "Sib2"), c(3,4))
#' x = simpleSim(x, 200, 1:2) # 200 equifrequent SNPs
#'
#' # Estimate the IBD coefficients for the siblings
#' est1 = IBDestimate(x, ids=c(3,4))
#'
#' # Estimate should be the same if pedigree structure is unknown
#' xlist = list(branch(x, 3), branch(x, 4))
#' est2 = IBDestimate(xlist, ids=c(3,4))
#' stopifnot(identical(est1, est2))
#'
#' # If the pedigree has plot.labels, they can be used as IDs
#' est3 = IBDestimate(x, ids=c("Sib1", "Sib2"))
#' stopifnot(identical(est1, est3))
#'
#' }
#'
#' @importFrom maxLik maxLik
#' @export
IBDestimate = function(x, ids, markers=NULL, start = c(0.99,0.001), tol=1e-7) {

  if (!requireNamespace("maxLik", quietly = TRUE)) {
    cat("Required package `maxLik` is not installed.")
    return(NULL)
  }

  single_linkdat = is.linkdat(x)
  if(single_linkdat) x = list(x)
  stopifnot(is.linkdat.list(x))

  if(is.null(markers))
      markers = seq_len(x[[1]]$nMark)
  SNPs = all(vapply(x[[1]]$markerdata[markers], attr, 'nalleles', FUN.VALUE=1)==2)

  if(is.vector(ids) && length(ids)==2)
      ids = rbind(ids)
  if(is.matrix(ids))
      ids = lapply(seq_len(nrow(ids)), function(i) pedlistMembership(x, ids[i,]))
  test_ids = is.list(ids) && all(sapply(ids, nrow) == 2)

  if(!test_ids) stop("Wrong format for parameter 'ids'.")

  # Optimize the above function in the triangle
  constraints = list(ineqA=matrix(c(1,0,-1,0,1,-1),3,2), ineqB=c(0,0,1)) # probability triangle

  res = lapply(ids, function(ids_df) {
      A = IBDest_getAlleleData(x, ids_df, markers)
      # Remove markers with missing alleles
      if(any(miss <- apply(A, 2, anyNA)))
          A = A[, !miss]

      # Likelihood function
      a=A[1,]; b=A[2,]; cc=A[3,]; d=A[4,]; pa=A[5,]; pb=A[6,]; pc=A[7,]; pd=A[8,]
      loglik_FUN = function(k) sum(log(.IBDlikelihood(k,a,b,cc,d,pa,pb,pc,pd)))

      # Optimise
      ML = maxLik(loglik_FUN, start=start, constraints=constraints, tol=tol)
      est = ML$estimate
      data.frame(ID1 = ids_df$orig.id[1], ID2 = ids_df$orig.id[2], Name1 = ids_df$plot.labels[1], Name2 = ids_df$plot.labels[2], N = ncol(A), k0 = est[1], k1=1-sum(est), k2=est[2], stringsAsFactors=FALSE)
  })

  do.call(rbind, res)
}


pedlistMembership = function(x, ids) {
    num = suppressWarnings(all(!is.na(as.numeric(ids))))
    if(num)
        ped_match = vapply(x, function(xx) ids %in% xx$orig.ids, FUN.VALUE=logical(length(ids)))
    else
        ped_match = vapply(x, function(xx) ids %in% xx$plot.labels, FUN.VALUE=logical(length(ids)))

    # Convert to matrix if length(ids==1)
    ped_match = rbind(ped_match, deparse.level=0)

    if(any(nomatch <- rowSums(ped_match) == 0))
        stop(paste("IDs not found in pedigrees:", paste(ids[nomatch], collapse=", ")))
    if(any(multimatch <- rowSums(ped_match) >1 ))
        stop(paste("IDs matching multiple pedigrees:", paste(ids[multimatch], collapse=", ")))

    pednr = apply(ped_match, 1, which)
    if(num) {
        orig.ids = ids
        plot.labels = vapply(seq_along(ids), function(i) {
            pedi = x[[pednr[i]]]
            if(is.null(pedi$plot.labels))
                return(NA_character_)
            int.id = match(ids[i], x[[pednr[i]]]$orig.ids)
            pedi$plot.labels[int.id]
        }, FUN.VALUE="a")
    }
    else {
        plot.labels = ids
        orig.ids = vapply(seq_along(ids), function(i) {
            pedi = x[[pednr[i]]]
            int.id = match(ids[i], x[[pednr[i]]]$plot.labels)
            pedi$orig.ids[int.id]
        }, FUN.VALUE=1)
    }
    data.frame(pednr=pednr, orig.ids=as.integer(orig.ids), plot.labels=plot.labels, stringsAsFactors=FALSE)
}

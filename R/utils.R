
.is.natural <- function(x) length(x) == 1 && is.numeric(x) && x == as.integer(x) && x > 0

.is.natural0 <- function(x) length(x) == 1 && is.numeric(x) && x == as.integer(x) && x >= 0

#' @export
on_failure(.is.natural) <- function(call, env) {
    paste0(deparse(call$x), " is not a positive integer")
}

#' @export
on_failure(.is.natural0) <- function(call, env) {
    paste0(deparse(call$x), " is not a non-negative integer")
}


.mysetdiff = function(x, y) unique.default(x[match(x, y, 0L) == 0L])
.myintersect = function(x, y) y[match(x, y, 0L)]

.internalID = function(x, orig.ids) {
    internal_ids = match(orig.ids, x$orig.ids)
    if (any(is.na(internal_ids))) 
        stop(paste("Indicated ID(s) not among original ID labels:", paste(orig.ids[is.na(internal_ids)], 
            collapse = ",")))
    internal_ids
}

.getSex = function(x, orig.ids) as.vector(x$pedigree[.internalID(x, orig.ids), "SEX"])

.comb2 = function(n) {
    if (n < 2) 
        return(matrix(nrow = 0, ncol = 2))
    v1 = rep.int(seq_len(n - 1), (n - 1):1)
    v2 = NULL
    for (i in 2:n) v2 = c(v2, i:n)
    cbind(v1, v2, deparse.level = 0)
}


.rand01 = function(n) sample.int(2, size = n, replace = T) - 1  #random 0/1 vector of length n.

.prettycat = function(v, andor) switch(min(len <- length(v), 3), toString(v), paste(v, collapse = " and "), 
    paste(paste(v[-len], collapse = ", "), andor, v[len]))



.generations = function(x) {
    # linkdat object
    max(vapply(unlist(.descentPaths(x, x$founders, original.ids = FALSE), recursive = F), length, 
        1))
}

.descentPaths = function(x, ids, original.ids = TRUE) {
    if (original.ids) 
        ids = .internalID(x, ids)
    offs = lapply(1:x$nInd, offspring, x = x, original.id = FALSE)
    lapply(ids, function(id) {
        res = list(id)
        while (TRUE) {
            newoffs = offs[vapply(res, function(path) path[length(path)], 1)]
            if (length(unlist(newoffs)) == 0) 
                break
            nextstep = lapply(1:length(res), function(r) if (length(newoffs[[r]]) == 0) 
                res[r] else lapply(newoffs[[r]], function(kid) c(res[[r]], kid)))
            res = unlist(nextstep, recursive = FALSE)
        }
        if (original.ids) 
            lapply(res, function(internal_vec) x$orig.ids[internal_vec]) else res
    })
}





# pedDistMatrix(x) computes a symmetrix integer matrix (d_ij), where d_ij is the pedigree
# distance (i.e., the length of the shortest pedigree path) between i and j. For example,
# the following relations have the indicated pedigree distances: parent/offspring 1,
# siblings 2, parents 2, uncle/niece 3, first cousins 4.  .pedDistMatrix <- function(x) { d
# = matrix(0, ncol=x$nInd, nrow=x$nInd) for (f in 1:x$nInd) { d[ f, offspring(x, f,
# original.id=FALSE) ] <- 1 d[ f, x$pedigree[f, c('FID', 'MID')] ] <- 1 } k=1 while
# (prod(d[col(d)!=row(d)])==0) { k=k+1 for (r in 1:nrow(cc <- .comb2(x$nInd))) { a=cc[r,1];
# b=cc[r,2] if (d[a,b] > 0) next if (any( d[a,]>0 & d[,b]>0 & d[a,]+d[,b] == k)) d[a,b] <-
# d[b,a] <- k } } d }

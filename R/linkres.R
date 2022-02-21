#' S3 methods for class 'linkres'.
#'
#' Functions for printing, summarizing and plotting the results of a linkage
#' analysis.
#'
#'
#' @param x,object a \code{linkres} object (normally produced by
#'   \code{\link{lod}} or \code{\link{merlin}}).
#' @param sort a logical, indicating if the data frame should be sorted
#'   according to map position.
#' @param chrom NULL, or a numeric containing chromosome numbers. In the latter
#'   case only results for the markers on the indicated chromosomes will be
#'   plotted.
#' @param ylim NULL, or a numeric of length 2: to be passed on to plot.default.
#' @param threshold a single numeric. A peak is defined as a regions of at least
#'   \code{width} consecutive markers LOD score above \code{threshold}.
#' @param width a single numeric.
#' @param physmap a matrix or data frame with three columns: Marker name,
#'   chromosome and physical position. This argument is optional.
#' @param \dots further arguments.
#' @seealso \code{\link{lod}}, \code{\link{merlin}}
#'
#' @examples
#'
#' x = linkdat(toyped, model=1)
#' lods = lod(x, theta='max')
#' summary(lods)
#' as.data.frame(lods)
#'
#' @name linkres
NULL

#' @rdname linkres
#' @export
print.linkres = function(x, ...) {
    rownames(x) = paste(rownames(x), ":", sep = "")
    if (attr(x, "analysis") == "mlink")
        rownames(x) = paste("theta=", rownames(x), sep = "")
    df = as.data.frame(unclass(x))
    print(df, print.gap = 2, ...)
}

#' @rdname linkres
#' @export
summary.linkres = function(object, ...) {
    x = object
    switch(attr(x, "analysis"), mlink = {
        mlods = apply(x, 2, function(i) if (all(is.na(i))) NA else max(i, 0, na.rm = TRUE))
        MX = max(mlods, na.rm = TRUE)
        cat("Max LOD score:", MX, "\n")
        cat("Achieved at marker(s):", names(mlods)[which(MX - mlods < 1e-04)], "\n")
    }, ilink = {
        lods = x["LOD", ]
        MX = max(lods, na.rm = TRUE)
        cat("Max LOD score:", MX, "\n")
        cat("Achieved at the following marker(s):\n")
        print(x[, which(MX - lods < 1e-04), drop = FALSE])
    })
}


#' @rdname linkres
#' @export
as.data.frame.linkres = function(x, ..., sort = TRUE) {
    map = attr(x, "map")
    map = map[map$MARKER %in% colnames(x), , drop = FALSE]
    if (sort)
        map = map[order(map$CHR, map$POS), , drop = FALSE]
    if (attr(x, "analysis") == "mlink")
        LOD = apply(x[, map$MARKER, drop = FALSE], 2, function(co) if (all(is.na(co)))
            NA else max(co, na.rm = TRUE)) else {
        LOD = x["LOD", ]
        x = x[2, , drop = FALSE]
    }
    lods = cbind(map, LOD = LOD, t(x[, map$MARKER, drop = FALSE]))
    lods
}



#' LOD score peaks
#'
#' Identify LOD score peaks
#'
#' The function first transforms \code{x} to a data frame (using
#' \code{\link{as.data.frame.linkres}} with \code{sort=T}. A peak is defined a
#' run of at least \code{width} consecutive markers with LOD score above or
#' equal to \code{threshold}. If possible, one flanking marker is included on
#' each side of the peak.
#'
#' @param x a \code{\link{linkres}} object
#' @param threshold a single numeric
#' @param width a positive integer
#' @return A list of data frames.
#' @seealso \code{\link{linkres}}, \code{\link{lod}}, \code{\link{merlin}},
#'
#' @examples
#'
#' ## minimal example
#' x = linkdat(toyped, model=1)
#' res = lod(x)
#' peak1 = lod.peaks(res, threshold=0)
#' peak2 = lod.peaks(res, threshold=0, width=2)
#' peak3 = lod.peaks(res, threshold=1)
#' stopifnot(length(peak1)==1, nrow(peak1[[1]])==1, length(peak2)==0, length(peak3)==0)
#'
#' @export
lod.peaks = function(x, threshold, width = 1) {
    # x et linkres objekt, eller data.frame med CHR, POS, LOD

    peak1chr = function(xchr, threshold, width) {
        # x must be sorted!
        rl = rle(xchr$LOD >= threshold)
        while (1) {
            short = rl$values & (rl$lengths < width)
            if (any(short)) {
                rl$values[short] <- FALSE
                rl = rle(inverse.rle(rl))
            } else break
        }
        xchr_nrow = nrow(xchr)
        if (!any(rl$values))
            return(NULL)
        start_ind = c(0, cumsum(rl$lengths))[which(rl$values)]
        stop_ind = start_ind + rl$lengths[rl$values] + 1  # plus 1 to compensate for endpoint[1]
        lapply(1:length(start_ind), function(i) {
            strt = start_ind[i]
            stp = stop_ind[i]
            telomeric = c(if (strt == 0) "start", if (stp > xchr_nrow) "end")
            if (length(telomeric) == 0)
                telomeric = "no"
            strt = max(strt, 1)
            stp = min(stp, xchr_nrow)
            structure(xchr[strt:stp, , drop = F], rownames = NULL, telomeric = telomeric)
        })
    }
    df = as.data.frame(x)
    df = df[!is.na(df$LOD), ]
    chrs = unique.default(df$CHR)
    res = list()
    for (chr in chrs) {
        dfchr = df[df$CHR == chr, ]
        res = c(res, peak1chr(dfchr, threshold = threshold, width = width))
    }
    res
}

#' @rdname linkres
#' @export
peakSummary = function(x, threshold, width = 1, physmap = NULL) {
    if (inherits(x, "linkres")) {
        if (is.null(threshold))
            stop("argument \"threshold\" is missing, with no default")
        x = lod.peaks(x, threshold, width)
    }
    dat = lapply(x, function(df) {
        n = nrow(df)
        chr = df$CHR[1]
        from_marker = df$MARKER[1]
        from_cm = df$POS[1]
        to_marker = df$MARKER[n]
        to_cm = df$POS[n]
        from_lod = df$LOD[1]
        max_lod = max(df$LOD)
        to_lod = df$LOD[n]
        L_cm = to_cm - from_cm
        telom = attr(df, "telomeric")
        data.frame(CHR = chr, FROM_MARKER = from_marker, TO_MARKER = to_marker, FROM_CM = from_cm,
            TO_CM = to_cm, TELOMERIC = telom, L_CM = L_cm, N = n, FROM_LOD = from_lod, MAX_LOD = max_lod,
            TO_LOD = to_lod, stringsAsFactors = F)
    })
    a = do.call(rbind, dat)
    if (!is.null(physmap)) {
        if (is.matrix(physmap))
            physmap = as.data.frame(physmap)
        if (is.character(physmap))
            physmap = read.table(physmap, header = T, as.is = T)
        if (is.data.frame(physmap)) {
            from_bp = physmap[match(a$FROM_MARKER, physmap[, 1]), 3]
            to_bp = physmap[match(a$TO_MARKER, physmap[, 1]), 3]
            a = cbind(a, FROM_BP = from_bp, TO_BP = to_bp)
            a = a[, c(1:5, 12:13, 6:11)]
        }
    }
    a
}


.getMap = function(x, markers = seq_len(x$nMark), na.action = 0, verbose = TRUE) {
    m = x$markerdata[markers]
    chrom = unlist(lapply(m, attr, "chrom"))
    marker = unlist(lapply(m, attr, "name"))
    pos = unlist(lapply(m, attr, "pos"))
    map = data.frame(CHR = chrom, MARKER = marker, POS = pos, stringsAsFactors = FALSE)
    if (na.action > 0) {
        nas = is.na(marker)
        map$MARKER[nas] = paste("M", markers[nas], sep = "")
    }
    if (na.action == 1) {
        nas2 = (is.na(chrom) | is.na(pos))
        if (all(nas2)) {
            if (verbose)
                cat("Warning: No map info given. Creating dummy map.\n")
            map$CHR = rep.int(1, x$nMark)
            map$POS = seq_len(x$nMark)
        }
    }
    # if(na.action ==2) nas2 = (is.na(chrom) | is.na(pos)) if(any(nas2)) { if(verbose)
    # cat('Warning: Deleting', sum(nas2), 'markers with missing map coordinates.\n') map =
    # map[!nas2, , drop=FALSE] }
    map
}

#' @rdname linkres
#' @export
plot.linkres = function(x, chrom = NULL, ylim = NULL, ...) {
    analysis = attr(x, "analysis")
    map = attr(x, "map")
    if (any(is.na(map$CHR))) {
        warning("Incomplete or missing map.")
        map = map[!is.na(map$CHR), ]
    }
    map = map[map$MARKER %in% colnames(x), , drop = FALSE]
    map = map[order(map$CHR, map$POS), , drop = FALSE]
    x = x[, match(map$MARKER, colnames(x)), drop = FALSE]

    if (!is.null(chrom)) {
        subindex = which(map$CHR %in% chrom)
        if (length(subindex) == 0)
            stop("No markers on indicated chromosome(s).")
        subx = structure(x[, subindex, drop = FALSE], analysis = analysis, map = map[subindex,
            , drop = F], class = "linkres")
        return(plot(subx, chrom = NULL, ylim = ylim, ...))
    }

    switch(analysis, mlink = {
        lds <- apply(x, 2, max)
        if (is.null(ylim)) ylim <- c(-1.2, max(c(3, lds), na.rm = T) + 0.3)
    }, ilink = {
        lds <- x["LOD", ]
        if (is.null(ylim)) ylim <- c(-0.5, max(c(3, lds), na.rm = T) + 0.3)
    }, )

    nM = ncol(x)
    pos = map$POS
    chr_br = which(map$CHR[-1] != map$CHR[-nM])
    for (b in chr_br) pos[(b + 1):nM] = pos[(b + 1):nM] + map$POS[b]  # NB: by now, map$POS != pos

    multichr = length(chr_br) > 0

    plot(pos, sapply(lds, max, ylim[1]), ylim = ylim, type = "l", lwd = 2, cex = 0.3, xlab = ifelse(multichr,
        "Chromosome", paste("Position (cM) on chromosome", map$CHR[1])), xaxt = ifelse(multichr,
        "n", par("xaxt")), ylab = "LOD score", ...)
    abline(h = 0, col = 2, lwd = 2)
    if (multichr)
        axis(1, at = c(0, pos[chr_br]), labels = map$CHR[c(chr_br, nM)], lwd.ticks = 2)
}

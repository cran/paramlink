#' Read dat file in LINKAGE format
#'
#' Converts dat files in LINKAGE format to dat/map/freq files in MERLIN format
#'
#' @param datfile character. The path to the dat file.
#' @param chrom integer chromosome number (needed to create the MERLIN map).
#' @param comment_string character indicating comments (which are removed before
#'   processing).
#' @param write_to a character prefix used for naming the output files, or NULL
#'   if no files should be written.
#' @return If \code{write_to} is NULL, a list of data.frames named \code{dat},
#'   \code{map} and \code{freq}.
#'
#' @examples
#'
#' # No example given.
#'
#' @export
readDatfile = function(datfile, chrom, comment_string="<<", write_to = NULL) {
    removePattern = paste0(comment_string, ".*")
    dat0 = sub(removePattern, "", readLines(datfile))
    dat = lapply(strsplit(dat0, split = " |\t"), function(v) v[v != ""])

    nMark = as.numeric(dat[[1]])[1] - 1
    xlinked = as.numeric(dat[[1]])[3]
    ordering = as.numeric(dat[[3]])
    markernames = sapply(dat[6 + xlinked + (1:nMark) * 2], "[", 4)
    pos = cumsum(as.numeric(dat[[length(dat) - 1]]))
    if (ordering[1] == 2)
        ordering = c(1, ordering)
    if (length(markernames) - length(pos) == 1)
        pos = c(0, pos)
    stopifnot(all(ordering == seq_len(nMark + 1)), length(dat) == 7 + xlinked + 2 * nMark +
        3)

    equal = (pos[-1] == pos[-length(pos)])
    k = 0
    for (i in 2:length(pos)) if (equal[i - 1])
        pos[i] = pos[i] + 1e-04 * (k <- k + 1) else k = 0  #if consecutive entries are equal, add 0.0001's.

    map = data.frame(CHR = chrom, MARKER = markernames, POS = pos, stringsAsFactors = F)

    freqlist = lapply(dat[7 + xlinked + (1:nMark) * 2], function(r) as.numeric(r))
    nalls = lengths(freqlist, use.names = F)
    L = sum(nalls) + length(nalls)
    cum = cumsum(c(1, nalls + 1))
    length(cum) = length(nalls)  #remove last
    col1 = rep("A", L)
    col1[cum] = "M"

    col2 = character(L)
    col2[cum] = markernames
    allalleles = unlist(lapply(nalls, seq_len))
    col2[-cum] = allalleles

    col3 = character(L)
    allfreqs = unlist(freqlist)
    col3[-cum] = format(allfreqs, scientifit = F, digits = 6)

    freq = cbind(col1, col2, col3, deparse.level = 0)

    merlindat = cbind(c("A", rep("M", nMark)), c("my_disease", markernames))
    if (!is.null(write_to)) {
        write.table(map, file = paste(write_to, "map", sep = "."), row.names = F, col.names = F,
            quote = F)
        write.table(merlindat, file = paste(write_to, "dat", sep = "."), row.names = F, col.names = F,
            quote = F)
        write.table(freq, file = paste(write_to, "freq", sep = "."), row.names = F, col.names = F,
            quote = F)
    }
    invisible(list(dat = merlindat, map = map, freq = freq))
}

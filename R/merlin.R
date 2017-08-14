#' MERLIN wrappers
#' 
#' Wrappers for the MERLIN software, providing multipoint LOD scores and other
#' computations on pedigrees with marker data. These functions require MERLIN 
#' to be installed and correctly pointed to in the PATH environment variable. 
#'
#' For these functions to work, MERLIN must be installed and the path to
#' merlin.exe included in the PATH variable. The \code{merlin} function is
#' first and foremost a wrapper to the parametric linkage functionality of
#' MERLIN. 
#'
#' By default the following MERLIN command is run (via a call to
#' \code{\link{system}}) after creating appropriate files in the current
#' working directory:
#' 
#' \preformatted{%
#' merlin -p _merlin.ped -d _merlin.dat -m _merlin.map -f _merlin.freq
#'        --model _merlin.model --tabulate --markerNames --quiet}
#' 
#' The resulting multipoint LOD scores are extracted from the output and
#' returned in R as a \code{\link{linkres}} object.
#' 
#' Additional command parameters can be passed on using the \code{options}
#' argument (this is simply pasted onto the MERLIN command, so dashes must be
#' included). For example, to obtain singlepoint LOD scores instead of
#' multipoint, set \code{options='--singlepoint'}. (The singlepoint scores
#' should agree with the results of \code{lod(x)}, except in cases where some
#' individuals have partial genotypes (see Examples).)
#' 
#' If \code{model=FALSE} the \code{--model merlin.model} part is removed from
#' the MERLIN command above. This is necessary for some calculations, e.g.
#' likelihoods (see Examples).
#' 
#' The \code{merlinUnlikely} function is a wrapper for MERLIN's '--error'
#' command. The syntax is similar to that of \code{\link{mendelianCheck}}.
#' 
#' @param x a \code{\link{linkdat}} object
#' @param markers an integer vector indicating which markers to use (default:
#' all).
#' @param model a logical: If TRUE (and x$model is not NULL), the file
#' 'merlin.model' is created and '--model merlin.model' is included to the
#' MERLIN command.
#' @param theta a numeric with values between 0 and 0.5: The recombination
#' value(s) for which the LOD score is computed.  The values of \code{theta}
#' are converted to centiMorgan positions using the Haldane map function and
#' included in the MERLIN command using the \code{--position} parameter. Works
#' only for single markers (i.e. \code{markers} must consist of a single
#' integer).
#' @param options a character with additional options to the MERLIN command.
#' See details.
#' @param verbose a logical: Show MERLIN output and other information, or not.
#' @param generate.files a logical. If TRUE, the files 'merlin.ped',
#' 'merlin.dat', 'merlin.map', 'merlin.freq' and (if \code{model=TRUE})
#' 'merlin.model' are created in the working directory.
#' @param cleanup a logical: Should the MERLIN files be deleted automatically?
#' @param logfile a character. If this is given, the MERLIN screen output will
#' be written to a file with this name.
#' @param remove a logical. If FALSE, the function returns the indices of
#' markers found to unlikely.  If TRUE, a new \code{linkdat} object is
#' returned, where the unlikely markers have been deleted.
#' @return If \code{model=TRUE}, a \code{\link{linkres}} object. Otherwise a
#' character containing the complete MERLIN output.
#' 
#' For \code{merlinUnlikely}, a numeric containing the indices of the unlikely,
#' or (if \code{remove=TRUE}) a new \code{linkdat} object where the unlikely
#' markers are removed.
#' @author Magnus Dehli Vigeland
#' @references \url{http://csg.sph.umich.edu/abecasis/Merlin/}
#' 
#' @examples
#' 
#' \dontrun{
#' x = linkdat(toyped, model=1)
#' x
#' 
#' # MERLIN treats partial genotypes (i.e. one known and one unknown allele) as missing:
#' lod_merlin = merlin(x)
#' lod_partial = lod(x)
#' x = modifyMarker(x, marker=1, ids=1, genotype=0)
#' lod_missing = lod(x)
#' stopifnot(lod_merlin == round(lod_missing, 4))
#' 
#' # Likelihood computation by MERLIN:
#' merlin(x, model=F, options='--lik')
#' }
#' 
#' @export
merlin = function(x, markers = seq_len(x$nMark), model = TRUE, theta = NULL, options = "", 
    verbose = FALSE, generate.files = TRUE, cleanup = generate.files, logfile = "") {
    
    clean = function(cleanup, verbose, files) if (cleanup) {
        unlink(files)
        if (verbose) 
            cat("Files successfully removed\n")
    }
    
    if (x$nMark == 0) 
        stop("No markers exist for this linkdat object.")
    if (model && is.null(x$model)) 
        stop("No model is set for this object")
    x = removeMarkers(x, seq_len(x$nMark)[-markers])
    map = .getMap(x, na.action = 1, verbose = F)
    mNames = map$MARKER
    
    extensions = c("ped", "dat", "map", "freq", if (model) "model")
    
    if (generate.files) {
        files = write.linkdat(x, prefix = "_merlin", what = extensions, merlin = TRUE)
        if (verbose) {
            cat("Files successfully generated:\n")
            print(files)
        }
    }
    
    options = paste(options, "--markerNames --quiet ")
    
    if (nonz_theta <- any(theta > 0)) {
        if (length(markers) > 1) {
            clean(cleanup, verbose, files)
            stop("Nonzero 'theta' values are possible only with a single marker.")
        }
        if (any(theta > 0.5)) {
            clean(cleanup, verbose, files)
            stop("Recombination fractions cannot exceed 0.5.")
        }
        pos = as.numeric(map[1, 3]) - 50 * log(1 - 2 * theta)  #Haldane's map: Converting rec.fractions to cM positions.
        options = paste(options, " --positions:", paste(pos, collapse = ","), sep = "")
    }
    
    program = if (identical(x$model$chrom, "X")) 
        "minx" else "merlin"
    command = paste(program, " -p _merlin.ped -d _merlin.dat -m _merlin.map -f _merlin.freq ", 
        if (model) 
            "--model _merlin.model --tabulate ", options, sep = "")
    
    if (verbose) 
        cat("\nExecuting the following command:\n", command, "\n\n", sep = "")
    
    merlinout = suppressWarnings(system(command, intern = T))
    clean(cleanup, verbose, files)
    if (nzchar(logfile)) 
        write(merlinout, logfile)
    if (any(substr(merlinout, 1, 11) == "FATAL ERROR")) {
        cat("\n====================================\n", paste(merlinout[-(2:10)], collapse = "\n"), 
            "====================================\n\n")
        return(invisible())
    }
    if (verbose) {
        cat("Merlin run completed\n")
        print(merlinout)
    }
    if (!is.na(skipped <- which(substr(merlinout, 3, 9) == "SKIPPED")[1])) 
        stop(paste(merlinout[c(skipped - 1, skipped)], collapse = "\n"))
    
    if (!model) 
        return(merlinout)
    
    ## Extract LOD scores
    res = read.table("merlin-parametric.tbl", sep = "\t", header = T, colClasses = c("numeric", 
        "numeric", "character", "NULL", "numeric", "NULL", "NULL"))  # chrom, pos, marker names and LOD.
    if (cleanup) 
        unlink("merlin-parametric.tbl")
    
    if (nonz_theta) {
        mlodsdim = c(length(theta), 1)
        dimnam = list(theta, mNames)
    } else {
        markernames = res$LABEL
        if (!all(markernames %in% mNames)) {
            markernames = paste(res$CHR, res$POS * 100, sep = "_")  # create new map with markernames formed as 'chr_pos'. Merlin weirdness: POS is in Morgans??
            map = data.frame(CHR = res$CHR, MARKER = markernames, POS = res$POS * 100)
        }
        mlodsdim = c(1, nrow(res))
        dimnam = list(0, markernames)
    }
    
    lodres = structure(res$LOD, dim = mlodsdim, dimnames = dimnam, analysis = "mlink", map = map, 
        class = "linkres")
    return(lodres)
}

#' @rdname merlin
#' @export
merlinUnlikely <- function(x, remove = FALSE, verbose = !remove) {
    merlin(x, model = F, options = "--error --prefix _merlinerror", verbose = verbose)
    err = read.table("_merlinerror.err", header = T, as.is = T)
    unlink("_merlinerror.err")
    if (verbose) 
        print(err)
    if (remove) 
        return(removeMarkers(x, markernames = err$MARKER)) else return(getMarkers(x, markernames = err$MARKER))
}

.merlin.unlikely = merlinUnlikely


.readMap = function(map, dat, freq, verbose, numerical = FALSE) {
    # TODO: numerical not used
    stopifnot(!is.null(map), !is.null(dat))
    if (is.character(map) && length(map) == 1) {
        rawmap = read.table(map, header = FALSE, comment.char = "", colClasses = "character")
        # If no number occurs in first entry, first row is assumed to be header.
        if (!any(1:9 %in% strsplit(rawmap[[1]][1], "")[[1]])) 
            rawmap = rawmap[-1, ]
    } else rawmap = map
    rawmap[[1]][rawmap[[1]] == "X"] = 23
    rawmap[[1]][rawmap[[1]] == "Y"] = 24
    rawmap[[1]][!rawmap[[1]] %in% 1:24] = NA
    mapnames = as.character(rawmap[[2]])
    map1 = matrix(as.numeric(c(rawmap[[1]], rawmap[[3]])), ncol = 2, dimnames = list(mapnames, 
        c("CHR", "POS")))
    
    if (is.character(dat) && length(dat) == 1) 
        rawdat = read.table(dat, header = FALSE, comment.char = "", colClasses = "character") else rawdat = dat
    datnames = as.character(rawdat[rawdat[[1]] == "M", 2])  #names of all markers in map
    
    if (is.character(freq) && length(freq) == 1) {
        rawfreq = as.matrix(read.table(freq, header = FALSE, colClasses = "character", fill = T))
        NROW = nrow(rawfreq)
        if (rawfreq[2, 1] == "F") {
            freqnames = rawfreq[seq(1, NROW, by = 2), 2]
            freqmatr = matrix(as.numeric(rawfreq[seq(2, NROW, by = 2), -1]), nrow = NROW/2, 
                dimnames = list(freqnames, NULL))
            isna = is.na(freqmatr)
            nalls = rowSums(!isna)
            allelmatr = col(freqmatr)
            allelmatr[isna] = NA
        } else {
            # if rawfreq[2,1]=='A' ... i.e. long format! See MERLIN tutorial about input files.
            Mrow = which(rawfreq[, 1] == "M")
            freqnames = rawfreq[Mrow, 2]
            
            nalls = c(Mrow[-1], NROW + 1) - Mrow - 1
            nM = length(Mrow)
            
            freqmatr = matrix(integer(), ncol = max(nalls), nrow = nM, dimnames = list(freqnames, 
                NULL))
            indexmatr = cbind(rep(seq_len(nM), times = nalls), unlist(lapply(nalls, seq_len)))
            freqmatr[indexmatr] = as.numeric(rawfreq[-Mrow, 3])
            
            allalleles = rawfreq[-Mrow, 2, drop = FALSE]
            numerical = !any(is.na(suppressWarnings(all_num <- as.numeric(allalleles))))
            if (numerical) {
                allelmatr = matrix(numeric(), ncol = max(nalls), nrow = nM, dimnames = list(freqnames, 
                  NULL))
                allelmatr[indexmatr] = all_num
            } else {
                allelmatr = matrix(character(), ncol = max(nalls), nrow = nM, dimnames = list(freqnames, 
                  NULL))
                allelmatr[indexmatr] = allalleles
            }
        }
        seqlist = lapply(1:max(nalls), seq_len)  # precomputing vectors to speed up mapinfo further down
    } else if (is.list(freq)) {
        stop("Unexpected frequency object in readmap()")
        freqlist = freq  #old
    } else freqmatr = NULL
    
    mapMatch = match(datnames, mapnames, nomatch = 0)
    if (verbose && any(NAs <- mapMatch == 0)) 
        cat("Deleting the following marker(s), which are not found in the map file:\n", paste(datnames[NAs], 
            collapse = "\n"), "\n", sep = "")
    
    if (is.null(freqmatr)) 
        annotations = lapply(seq_along(datnames), function(i) {
            if ((mapmatch = mapMatch[i]) == 0) 
                return(NULL)
            list(chrom = map1[mapmatch, 1], pos = map1[mapmatch, 2], name = datnames[i])
        }) else {
        freqMatch = match(datnames, freqnames, nomatch = 0)
        annotations = lapply(seq_along(datnames), function(i) {
            if ((mapmatch <- mapMatch[i]) == 0) 
                return(NULL)
            if ((fmatch <- freqMatch[i]) == 0) 
                alleles = afreq = NULL else {
                sq = seqlist[[nalls[fmatch]]]
                alleles = allelmatr[fmatch, sq]  ###TODO: CHECK THIS!!
                afreq = as.vector(freqmatr[fmatch, sq])
            }
            list(alleles = alleles, afreq = afreq, chrom = map1[mapmatch, 1], pos = map1[mapmatch, 
                2], name = datnames[i])
        })
    }
    annotations
}

#' Linkdat objects
#'
#' Functions to create and display 'linkdat' objects.
#'
#' The file (or matrix or data.frame) \code{ped} must describe one or several
#' pedigrees in standard LINKAGE format, i.e. with the following columns in
#' correct order:
#'
#' 1 Family id (optional) (FAMID)
#'
#' 2 Individual id (ID),
#'
#' 3 Father id (FID),
#'
#' 4 Mother id (MID),
#'
#' 5 Gender (SEX): 1 = male, 2 = female,
#'
#' 6 Affection status (AFF): 1 = unaffected, 2 = affected, 0 = unknown,
#'
#' 7 First allele of first marker,
#'
#' 8 Second allele of first marker,
#'
#' 9 First allele of second marker,
#'
#' a.s.o.
#'
#' Only columns 2-6 are mandatory. The first column is automatically interpreted
#' as family id if it has repeated elements.
#'
#' Internally the individuals are relabeled as 1,2,..., but this should rarely
#' be of concern to the end user. Some pedigree checking is done, but it is
#' recommended to plot the pedigree before doing any analysis.
#'
#' Details on the formats of map, dat and frequency files can be found in the
#' online MERLIN tutorial: \url{http://csg.sph.umich.edu/abecasis/Merlin/}
#'
#' A singleton is a special \code{linkdat} object whose pedigree contains 1
#' individual. The class attribute of a singleton is \code{c('singleton',
#' 'linkdat')}
#'
#' @param ped a matrix, data.frame or a character with the path to a pedigree
#'   file in standard LINKAGE format. (See details)
#' @param model either a \code{linkdat.model} object (typically \code{y$model}
#'   for some linkdat object \code{y}), or a single integer with the following
#'   meaning: 1 = autosomal dominant; 2 = autosomal recessive; 3 = X-linked
#'   dominant; 4 = X-linked recessive. In each of these cases, the disease is
#'   assumed fully penetrant and the disease allele frequency is set to 0.00001.
#'   If \code{model=NULL}, no model is set.
#' @param map a character with the path to a map file in MERLIN format, or NULL.
#'   If non-NULL, a dat file must also be given (next item).
#' @param dat a character with the path to a dat file in MERLIN format, or NULL.
#'   (Only needed if \code{map} is non-NULL.)
#' @param freq a character with the path to a allele frequency file in MERLIN
#'   (short) format, or NULL. If NULL, all markers are interpreted as
#'   equifrequent.
#' @param annotations a list (of the same length and in the same order as the
#'   marker columns in \code{x}) of marker annotations. If this is non-NULL,
#'   then all of \code{map, dat, freq} should be NULL.
#' @param missing the character (of length 1) used for missing alleles. Defaults
#'   to '0'.
#' @param header a logical, relevant only if \code{ped} points to a ped file: If
#'   TRUE, the first line of the ped file is skipped.
#' @param checkped a logical. If FALSE, no checks for pedigree errors are
#'   performed.
#' @param verbose a logical: verbose output or not.
#' @param id,sex single numerics describing the individual ID and gender of the
#'   singleton.
#' @param markers a numeric indicating which markers should be included/printed.
#' @param x,object a \code{linkdat} object.
#' @param famid a numeric: the family ID of the singleton.
#' @param subset a numeric containing the individuals in the sub-pedigree to be
#'   extracted. NB: No pedigree checking is done here, so make sure the subset
#'   form a meaningful, closed pedigree.
#' @param prefix a character string giving the prefix of the files. For
#'   instance, if \code{prefix='fam1'} and \code{what=c('ped', 'map')}, the
#'   files 'fam1.ped' and 'fam1.map' will be created.
#' @param what a character vector forming a subset of c('ped', 'map', 'dat',
#'   'freq', 'model'), indicating which files should be created. All files are
#'   written in MERLIN style (but see the next item!)
#' @param merlin a logical. If TRUE, the marker alleles are relabeled to
#'   1,2,..., making sure that the generated files are readable by MERLIN (which
#'   does not accept non-numerical allele labels in the frequency file.) If
#'   FALSE (the default) the allele labels are unchanged. In this case, \code{x}
#'   should be exactly reproducible from the files. (See examples.)
#' @param \dots further arguments.
#' @return A \code{linkdat} object, or a list of \code{linkdat} objects. A
#'   linkdat object is essentially a list with the following entries, some of
#'   which can be NULL.  \item{pedigree }{\code{data.frame} with 5 columns (ID,
#'   FID, MID, SEX, AFF) describing the pedigree in linkage format. (NB:
#'   Internal labeling used.)} \item{orig.ids}{the original individual id
#'   labels.} \item{nInd}{the number of individuals in the pedigree.}
#'   \item{founders}{vector of the founder individuals. (NB: Internal labeling
#'   used.)} \item{nonfounders}{vector of the nonfounder individuals (NB:
#'   Internal labeling used.)} \item{hasLoops}{a logical: TRUE if the pedigree
#'   is inbred.} \item{subnucs}{list containing all (maximal) nuclear families
#'   in the pedigree. Each nuclear family is given as a vector of the form
#'   c(pivot, father, mother, child1, ...), where the pivot is either the id of
#'   the individual linking the nuclear family to the rest of the pedigree, or 0
#'   if there are none. (NB: Internal labeling used.)} \item{markerdata}{a list
#'   of \code{\link{marker}} objects.} \item{nMark}{the number of markers.}
#'   \item{available}{a numeric vector containing IDs of available individuals.
#'   Used for simulations and plots.} \item{model}{a \code{linkdat.model}
#'   object, essentially a list containing the model parameters. See
#'   \code{\link{setModel}} for details.} \item{loop_breakers}{a matrix with
#'   original loop breaker ID's in the first column and their duplicates in the
#'   second column. This is set by \code{\link{breakLoops}}.}
#' @seealso \code{\link{pedCreate}}, \code{\link{pedModify}},
#'   \code{\link{pedParts}}, \code{\link{setModel}}
#'
#' @examples
#'
#' x = linkdat(toyped, model=1)
#' x
#' summary(x)
#'
#' #### test read/write:
#' x = modifyMarker(x, 1, alleles=c('B','C'), afreq=c(.9, .1), chrom=2, name='SNP1', pos=123)
#' write.linkdat(x, prefix='toy')
#' y = linkdat('toy.ped', map='toy.map', dat='toy.dat', freq='toy.freq', model=1)
#' unlink(c('toy.ped', 'toy.map', 'toy.dat', 'toy.freq', 'toy.model'))
#' stopifnot(isTRUE(all.equal(x,y)))
#'
#' #### test singletons:
#' w = singleton(id=3, sex=2)
#' T1 = all.equal(w, linkdat(ped=rbind(c(3,0,0,2,1))))
#' w = markerSim(w, N=5, alleles=2, afreq=c(0.1,.9))
#' T2 = all.equal(w, relabel(relabel(w, 10), 3))
#' T3 = all.equal(w, swapSex(swapSex(w, 3), 3))
#' T4 = all.equal(w, swapAff(swapAff(w, 3), 3))
#' stopifnot(T1, T2, T3, T4)
#'
#' #### several ways of creating the same linkdat object:
#' alleles = c(157,160,163)
#' afreq = c(0.3, 0.3, 0.4)
#' gt10 = c(160, 160)
#' gt14 = c(160, 163)
#'
#' z1 = relabel(addOffspring(nuclearPed(1), father=3, noffs=1, aff=2), 10:14)
#' z1 = addMarker(z1, marker(z1, 10, gt10, 14, gt14, alleles=alleles, afreq=afreq))
#' z1 = setModel(z1, 2)
#'
#' z2 = addParents(relabel(nuclearPed(1), 12:14), 12, father=10, mother=11)
#' z2 = addMarker(z2, rbind(gt10, 0, 0, 0, gt14), alleles=alleles, afreq=afreq)
#' z2 = setModel(swapAff(z2, 14), 2)
#'
#' z3 = linkdat(data.frame(ID=10:14, FID=c(0,0,10,0,12), MID=c(0,0,11,0,13),
#'              SEX=c(1,2,1,2,1), AFF=c(1,1,1,1,2),
#'              M=c('160/160', '0/0', '0/0', '0/0', '160/163')), model=2)
#' z3 = modifyMarker(z3, 1, alleles=alleles, afreq=afreq)
#'
#' write.linkdat(z1, prefix='test')
#' z4 = linkdat('test.ped', map='test.map', dat='test.dat', freq='test.freq',
#'              model=2)
#' z4 = modifyMarker(z4, 1, alleles=alleles, chrom=NA, pos=NA, name=NA)
#'
#' write.linkdat(z1, prefix='test', merlin=TRUE)
#' z5 = linkdat('test.ped', map='test.map', dat='test.dat', freq='test.freq',
#'              model=2)
#' z5 = modifyMarker(z5, 1, alleles=alleles, chrom=NA, pos=NA, name=NA)
#'
#' stopifnot(isTRUE(all.equal(z1,z2)), isTRUE(all.equal(z1,z3)),
#'           isTRUE(all.equal(z1,z4)), isTRUE(all.equal(z1,z5)))
#' unlink(c('test.ped', 'test.map', 'test.dat', 'test.freq', 'test.model'))
#'
#' @export
linkdat = function(ped, model = NULL, map = NULL, dat = NULL, freq = NULL, annotations = NULL,
    missing = 0, header = FALSE, checkped = TRUE, verbose = TRUE, ...) {

    subnucs <- function(ped) {
        # output: peeling order of nuclear subfamilies. Format for each nuc:
        # c(pivot,father,mother,offsp1,..), where pivot=0 for the last nuc.
        if (nrow(ped) == 1)
            return(list())
        parents = unique(ped[, 2:3])
        parents = parents[-match(0, parents[, 1]), , drop = FALSE]
        list1 = lapply(nrow(parents):1, function(i) {
            par = parents[i, ]
            list(father = par[[1]], mother = par[[2]], offspring = as.vector(ped[, 1])[which(ped[,
                2] == par[[1]] & ped[, 3] == par[[2]], useNames = FALSE)])
        })  #listing all nucs
        res = list()
        i = 1
        k = 1
        while (length(list1) > 1) {
            if (i > length(list1))
                return(FALSE)
            sub = list1[[i]]
            subvec = unlist(sub)
            links = subvec[subvec %in% unlist(list1[-i])]
            if (length(links) == 1) {
                res[[k]] <- c(sub, list(pivot = as.numeric(links), pivtype = match(links, c(sub[["father"]],
                  sub[["mother"]]), nomatch = 3)))
                list1 <- list1[-i]
                k <- k + 1
                i <- 1
            } else i <- i + 1
        }
        res[[k]] <- c(list1[[1]], list(pivot = 0, pivtype = 0))  #final nuclear
        res
    }

    if (is.linkdat(ped))
        return(ped)
    if (is.character(ped) && length(ped) == 1) {
        skip = as.integer(header)
        first = scan(ped, what = "", skip = skip, nlines = 1, quiet = TRUE, ...)
        ncols = length(first)
        .numerical = !any(is.na(suppressWarnings(as.numeric(first))))
        ped = scan(ped, what = "", skip = skip, quiet = TRUE, ...)
        ped = matrix(ped, ncol = ncols, byrow = TRUE)
    }
    ped = as.matrix(ped)  # if ped is a data frame
    if (!exists(".numerical"))
        .numerical = !is.numeric(ped)  # looks wrong, but makes sense in next line: If ped is a character matrix
    numerical = .numerical && !any(is.na(suppressWarnings(ped_num <- as.numeric(ped))))
    if (numerical) {
        dim(ped_num) = dim(ped)
        ped = ped_num
    }

    nrows = nrow(ped)
    if (nrows == 0)
        stop("Empty pedigree.")

    if (!is.null(map)) {
        if (is.null(dat))
            stop("The 'map' and 'dat' arguments must either both be NULL, or both non-NULL")
        annotations = .readMap(map, dat, freq, verbose, numerical)
    } else if (is.character(annotations))
        annotations = list(alleles = switch(annotations, snp12 = 1:2, snpAB = c("A", "B")),
            afreq = c(0.5, 0.5))

    if (length(unique(ped[, 1])) < nrows | all(gsub(" ", "", ped[, 2], fixed = TRUE) != 0)) {
        # added second test, in case all rows are singletons
        famids = unique.default(ped[, 1])  # these are characters here
        if (length(famids) > 1)
            return(lapply(famids[order(as.numeric(famids))], function(fam) linkdat(ped[ped[,
                1] == fam, , drop = F], model = model, annotations = annotations, verbose = verbose,
                missing = missing))) else {
            famid = as.numeric(ped[1, 1])
            ped = ped[, -1, drop = F]
        }
    } else if (nrows == 1 && all(ped[, 3:4] == 0)) {
        # singleton with famid!
        famid = as.numeric(ped[1, 1])
        ped = ped[, -1, drop = F]
    } else famid = 1
    if (verbose)
        cat("Family ID: ", famid, ".\n", sep = "")

    if (ncol(ped) < 5)
        stop("Too few columns: ID, FID, MID, SEX and AFF are mandatory.")
    pedcols = ped[, 1:5]
    if (!numerical) {
        if (any(grepl("[^0-9 ]", pedcols)))
            stop("Pedigree columns must be numeric.")
        pedcols = as.numeric(pedcols)
    }
    pedigree = matrix(pedcols, ncol = 5, dimnames = list(NULL, c("ID", "FID", "MID", "SEX",
        "AFF")))

    if (checkped)
        .checkped(pedigree)

    orig.ids = as.vector(pedigree[, 1])
    nInd = nrows
    pedigree = relabel(pedigree, new = 1:nInd)
    if (verbose) {
        if (nInd == 1)
            cat(sprintf("Singleton %s, individual ID = %d.\n", ifelse(pedigree[, "SEX"] ==
                1, "male", "female"), orig.ids)) else cat(nInd, "individuals.\n")
        affs = sum(pedigree[, "AFF"] == 2)
        if (affs > 0)
            cat(sprintf("%d affected, %d non-affected.\n", affs, nInd - affs))
    }
    founders = as.integer(which(pedigree[, "FID"] == 0))
    nonfounders = as.integer(which(pedigree[, "FID"] > 0))

    #---peeling order of nuclear subfamilies---
    if (nInd > 1) {
        subnucs = subnucs(pedigree)
        hasLoops = is.logical(subnucs) && !subnucs
        if (verbose)
            if (hasLoops)
                cat("Loop(s) detected.\n") else if (is.list(subnucs))
                cat(sprintf("%d nuclear %s.\n", length(subnucs), ifelse(length(subnucs) ==
                  1, "subfamily", "subfamilies")))
        class = "linkdat"
    } else {
        class = c("singleton", "linkdat")
        hasLoops = FALSE
        subnucs = NULL
    }

    #---creation of linkdat object---
    obj = structure(list(pedigree = pedigree, famid = famid, orig.ids = orig.ids, nInd = nInd,
        founders = founders, nonfounders = nonfounders, hasLoops = hasLoops, subnucs = subnucs),
        class = class)

    #---adding markers---
    obj = setMarkers(obj, ped[, -(1:5), drop = F], annotations = annotations, missing = missing)
    if (verbose)
        if (obj$nMark == 1)
            cat("1 marker.\n") else cat(obj$nMark, "markers.\n")

    #----adding model----
    if (!is.null(model))
        obj = setModel(obj, model = model)

    if (verbose)
        cat("\n")
    invisible(obj)
}


#' @export
#' @rdname linkdat
singleton = function(id, sex = 1, famid = 1, verbose = FALSE, ...) linkdat(ped = rbind(c(famid,
    id, 0, 0, sex, 1)), verbose = verbose, ...)


.checkped = function(p) {
    # p a numeric matrix with 5 columns
    ID = p[, "ID"]
    FID = p[, "FID"]
    MID = p[, "MID"]
    SEX = p[, "SEX"]
    AFF = p[, "AFF"]

    # singletons:
    if (nrow(p) == 1) {
        if (FID != 0 || MID != 0)
            stop("Singleton error: FID and MID must both be zero.")
        if (!SEX %in% 1:2)
            stop("Singleton error: Unknown sex.")
        return()
    }
    # real pedigrees:
    if (all(c(FID, MID) == 0))
        stop("Pedigree is not connected.")

    fatherErr = !FID %in% c(0, ID)
    motherErr = !MID %in% c(0, ID)
    self_ancest = sapply(seq_along(ID), function(i) ID[i] %in% ancestors(p, ID[i]))
    quick.check <- (all(SEX %in% 1:2) && all(AFF %in% 0:2) && all((FID > 0) == (MID > 0)) &&
        !any(duplicated(ID)) && !any(fatherErr) && !any(motherErr) && all(SEX[match(FID[FID !=
        0], ID)] == 1) && all(SEX[match(MID[MID != 0], ID)] == 2) && !any(self_ancest))

    if (quick.check)
        return()  #if all tests are passed


    for (i in seq_along(ID)) {
        if (!SEX[i] %in% 1:2)
            cat("Individual ", ID[i], ": SEX must be either 1 (male) or 2 (female).\n", sep = "")
        if (!AFF[i] %in% 0:2)
            cat("Individual ", ID[i], ": Affection status must be either 0 (unknown), 1 (non-affected) or 2 (affected).\n",
                sep = "")
        if ((FID[i] > 0) != (MID[i] > 0))
            cat("Individual ", ID[i], ": Only one parent in the pedigree is not allowed. Either both parents or none must be specified.\n",
                sep = "")
        if (i > 1 && ID[i] %in% ID[1:(i - 1)]) {
            cat("Individual ", ID[i], ": ID not unique.\n", sep = "")
            next
        }
        if (fatherErr[i])
            cat("Individual ", ID[i], ": Father's ID (", FID[i], ") does not appear in ID column.\n",
                sep = "") else if (FID[i] != 0 && SEX[match(FID[i], ID)] != 1)
            cat("Individual ", ID[i], ": Father is not male.\n", sep = "")
        if (motherErr[i])
            cat("Individual ", ID[i], ": Mother's ID (", MID[i], ") does not appear in ID column.\n",
                sep = "") else if (MID[i] != 0 && SEX[match(MID[i], ID)] != 2)
            cat("Individual ", ID[i], ": Mother is not female.\n", sep = "")
        if (self_ancest[i])
            cat("Individual ", ID[i], " is ", switch(SEX[i], "his", "her"), " own ancestor.\n",
                sep = "")
    }
    stop("Pedigree errors detected.")
}

#' @export
#' @rdname linkdat
print.linkdat = function(x, ..., markers) {
    if (missing(markers))
        marker.nos = seq_len(min(x$nMark, 5)) else {
        if (length(markers) > 0 && max(markers) > x$nMark)
            stop("Nonexisting marker(s) indicated")
        marker.nos = markers
    }
    datafr = as.data.frame(x, markers = marker.nos, sep = "/", missing = "-", singleCol = TRUE)
    print(datafr, ...)
    if (missing(markers) && x$nMark > 5)
        cat("\nOnly first 5 markers are shown. Use option 'markers=' to print specified markers.\n")
    invisible(datafr)
}

#' @export
#' @rdname linkdat
summary.linkdat = function(object, ...) {
    x <- object
    cat("Pedigree:\n---------\n")
    cat(x$nInd, "individuals\n")
    cat(length(x$founders), "founders,", length(x$nonfounders), "nonfounders; bit size =",
        2 * length(x$nonfounders) - length(x$founders), "\n")
    if ((ant <- length(x$subnucs)) > 0)
        cat(ant, "nuclear", ifelse(ant == 1, "subfamily", "subfamilies"), "\n")
    aff = x$pedigree[, "AFF"]
    if (all(aff == 1))
        cat("No pedigree members affected by disease\n") else cat(sum(aff == 2), "affected by disease,", sum(aff == 1), "unaffected,", sum(aff ==
        0), "with unknown affection status\n")

    cat("\nMarker data:\n------------\n", x$nMark, ifelse(x$nMark == 1, " marker ", " markers "),
        "in total\n", sep = "")
    if (x$nMark > 0) {
        miss = which(rowSums(m <- do.call(cbind, x$markerdata)) == 0)
        cat(length(miss), "individuals with no available genotypes")
        if (length(miss) > 0) {
            cat(":", paste(x$orig.ids[miss], collapse = ", "), "\n")
            cat(round(sum(m[-miss, ] == 0)/length(m[-miss, ]) * 100, 2), "% missing alleles (excluding ungenotyped individuals)\n")
        } else cat("\n", sum(m == 0)/length(m) * 100, "% missing alleles\n", sep = "")
        cat("\nChromosome distribution of markers:\n")
        chrtbl = table(sapply(x$markerdata, attr, "chrom"), useNA = "ifany")
        names(chrtbl)[is.na(names(chrtbl))] = "unknown"
        for (i in seq_along(chrtbl)) cat(" chromosome ", names(chrtbl)[i], ": ", chrtbl[i],
            ifelse(chrtbl[i] == 1, " marker\n", " markers\n"), sep = "")

        cat("\nAllele number distribution:\n")
        ntbl = table(unlist(lapply(x$markerdata, attr, "nalleles")))
        for (i in seq_along(ntbl)) cat(" ", names(ntbl)[i], " alleles", ": ", ntbl[i], ifelse(ntbl[i] ==
            1, " marker\n", " markers\n"), sep = "")
    }


    cat("\nModel parameters:\n-----------------\n")
    if (is.null(x$model))
        cat("No model parameters set\n") else print(x$model)
}

#' @export
#' @rdname linkdat
write.linkdat = function(x, prefix = "", what = c("ped", "map", "dat", "freq", "model"), merlin = FALSE) {
    generated.files = character(0)
    if (merlin) {
        pmatr = as.matrix(x)
        markerattr = attr(pmatr, "markerattr")
    }

    if (any(c("map", "dat", "freq") %in% what))
        .map = .getMap(x, na.action = 1, verbose = F)

    if ("ped" %in% what) {
        if (!merlin)
            pmatr = cbind(x$famid, relabel(x$pedigree, x$orig.ids),
                          .prettyMarkers(x$markerdata, singleCol = FALSE))
        write(t(pmatr), pedname <- paste(prefix, "ped", sep = "."), ncolumns = 6 + 2 * x$nMark)
        generated.files = c(generated.files, pedname)
    }

    if ("model" %in% what && !is.null(x$model)) {
        dfreq = format(x$model$dfreq, scientific = F, decimal.mark = ".")
        # if X-linked: using the female penetrances. This works in MERLIN/MINX as long as male and
        # female penetrances are equal.
        penetrances = if (x$model$chrom == "X")
            x$model$penetrances$female else x$model$penetrances
        penets = paste(format(penetrances, scientific = F, decimal.mark = "."), collapse = ",")
        .model = c("my_disease", dfreq, penets, "my_model")
        write(.model, modelname <- paste(prefix, "model", sep = "."), sep = " \t", ncolumns = 4)
        generated.files = c(generated.files, modelname)
    }

    if ("map" %in% what) {
        write.table(.map, mapname <- paste(prefix, "map", sep = "."), col.names = F, row.names = F,
            quote = F)
        generated.files = c(generated.files, mapname)
    }

    if ("dat" %in% what) {
        .dat = cbind(code = c("A", rep("M", nrow(.map))), value = c("my_disease", .map$MARKER))
        write.table(.dat, datname <- paste(prefix, "dat", sep = "."), col.names = F, row.names = F,
            quote = F)
        generated.files = c(generated.files, datname)
    }

    if ("freq" %in% what) {
        if (!merlin)
            markerattr = lapply(x$markerdata, attributes)
        nalls = unlist(lapply(markerattr, function(at) at$nalleles))
        L = sum(nalls) + length(nalls)

        cum = cumsum(c(1, nalls + 1))
        length(cum) = length(nalls)  #remove last

        col1 = rep("A", L)
        col1[cum] = "M"

        col2 = character(L)
        col2[cum] = .map$MARKER
        if (merlin)
            allalleles = unlist(lapply(nalls, seq_len)) else allalleles = unlist(lapply(markerattr, function(at) at$alleles))
        col2[-cum] = allalleles

        col3 = character(L)
        allfreqs = unlist(lapply(markerattr, function(at) at$afreq))
        col3[-cum] = format(allfreqs, scientifit = F, digits = 6)

        .freq = cbind(col1, col2, col3)
        write.table(.freq, freqname <- paste(prefix, "freq", sep = "."), col.names = F, row.names = F,
            quote = F)
        generated.files = c(generated.files, freqname)
    }

    invisible(generated.files)
}

#' @export
#' @rdname linkdat
subset.linkdat <- function(x, subset = x$orig.ids, ..., markers = seq_len(x$nMark)) {
    x = removeMarkers(x, setdiff(seq_len(x$nMark), markers))
    xframe = as.matrix(x)

    newfr = xframe[xframe[, "ID"] %in% subset, , drop = F]
    newfr[!(newfr[, "FID"] %in% subset), "FID"] = 0  # set FID=0 if father is not in subset
    newfr[!(newfr[, "MID"] %in% subset), "MID"] = 0  # set MID=0 if mother is not in subset

    restore_linkdat(newfr, attributes(xframe))
}

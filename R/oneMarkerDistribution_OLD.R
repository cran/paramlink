

# .oneMarkerDistribution_OLD <- function(x, ids, partialmarker, theta=0.5, loop_breakers=NULL, verbose=TRUE) {
	# if (theta < 0.5 && is.null(x$model)) stop("Model must be set first. Use setModel().") 
	# nids = length(ids)
	# if (is.numeric(partialmarker) && length(partialmarker)==1) {
		# if (partialmarker < 1 || partialmarker > x$nMark) stop("Indicated marker does not exist.")
		# else  partialmarker = x$markerdata[[partialmarker]]
	# }
	# if (class(partialmarker) != "marker") stop("The 'partialmarker' must be either a single integer or a 'marker' object.")
	# x = setMarkers(x, m <- partialmarker)
	# afreq = attr(m, 'afreq')
	# alleles = attr(m, "alleles")
	# if(is.null(x$model)) {
		# if(identical(23L, as.integer(attr(m, 'chrom')))) x = setModel(x, 4, penetrances=list(male=c(.5, .5), female=c(.5,.5,.5)))
		# else x = setModel(x, 1, penetrances=c(.5,.5,.5))
	# }
	# chrom = x$model$chrom
	# SEX = x$pedigree[,'SEX']
	
	# if(verbose) {
		# cat("Assuming that the marker is", ifelse(chrom=="AUTOSOMAL", "autosomal.", "X-linked."), "\n")
		# if(theta < 0.5)	cat("Recombination rate between marker and disease locus:", theta,".\n", sep="")
		# else if(any(x$pedigree[,'AFF'] != 1)) cat("Assuming theta=0.5, i.e that the marker is independent of disease status.\n")
		
		# cat("\nPartial marker data:\n")
		# print(data.frame(ID=x$orig.ids, GENO=.prettyMarkers(list(m), missing="-", singleCol=TRUE, chrom=chrom, sex=SEX)), row.names=FALSE)
		# cat("\nAllele frequencies:\n")
		# print(structure(afreq, names=alleles))
	# }
		
	# if (x$hasLoops)	{		
		# if(is.null(lb <- loop_breakers)) 	 stop("The pedigree has loops. Please indicate loop breakers.")
		# if(verbose) cat(ifelse(length(lb)==1, "\nBreaking loop at individual", "\nBreaking loops at individuals"), .prettycat(lb, "and"), "\n")
		# x = breakLoops(x, lb)
		# m = x$markerdata[[1]]
	# }
	
	# int.ids = .internalID(x, ids)
	# heterogenos = .comb2(nall <- attr(m, 'nalleles'))
	# allgenos = rbind(cbind(1:nall,1:nall), heterogenos)
	
	# switch(chrom,
	# AUTOSOMAL = {
		# genos = lapply(int.ids, function(id) { geno = m[id, ];  switch(sum(geno == 0) + 1, rbind(geno), cbind(geno[geno!=0], 1:nall), allgenos)	})
		# geno.names = lapply(genos, function(matr) paste(alleles[matr[,1]], alleles[matr[,2]], sep=""))
		# indexgrid = .my.grid(lapply(genos, function(matr) 1:nrow(matr))) 
		# probs = array(0, dim=sapply(geno.names, length), dimnames=geno.names)
		# for (r in 1:nrow(indexgrid)) {
			# indices = indexgrid[r,]
			# for(i in 1:nids) 
				# m[int.ids[i], ] = genos[[i]][indices[i], ]
			# #start.dat = .startdata_MD(x, m)
			# probs[rbind(indices)] = likelihood(x, locus1=m, locus2="disease", theta)  #.likelihood_multi(x, theta, dat=start.dat)
		# }
	# },
	# X = {
		# genos = lapply(int.ids, function(id) { geno = m[id, ];	switch(SEX[id], {if(geno[1]==0)  cbind(1:nall,1:nall) else rbind(geno, deparse.level=0)}, switch(sum(geno == 0) + 1, rbind(geno), cbind(geno[geno!=0], 1:nall), allgenos)) })
		# geno.names = lapply(seq_along(int.ids), function(i) { matr = genos[[i]]; switch(SEX[int.ids[i]], alleles[matr[,1]], paste(alleles[matr[,1]], alleles[matr[,2]], sep="")) })
		# indexgrid = .my.grid(lapply(genos, function(matr) 1:nrow(matr))) 
		# probs = array(0, dim=sapply(geno.names, length), dimnames=geno.names)
		# for (r in 1:nrow(indexgrid)) {
				# indices = indexgrid[r,]
				# for(i in 1:nids) 
					# m[int.ids[i], ] = genos[[i]][indices[i], ]
				# #start.dat = .startdata_MD_X(x, m)
				# probs[rbind(indices)] = likelihood(x, locus1=m, locus2="disease", theta) #.likelihood_multi_X(x, theta, dat=start.dat)
		# }
	# })

	# if (sum(probs)==0)
		# warning("\nAll genotype probabilities zero. Mendelian error?") 
	# res = probs/sum(probs)
	# if(verbose) {
		# cat(ifelse(nids==1, "\nGenotype probability distribution for individual ",
				# "\nJoint genotype probability distribution for individuals "), .prettycat(ids,"and"), ":\n", sep="")
		# print(round(res,4))
		# return(invisible(res))
	# }
	# else res
# }


# ##ANNET GAMMALT

# # .startdata_M_X = function(x, marker, eliminate=0) {
	# # marker = .reduce_alleles(marker)
	# # glist = .eliminate_X(x, .build_genolist_X(x, marker), attr(marker, 'nalleles'), repeats=eliminate)
	# # if (attr(glist, 'impossible')) {dat=list(); attr(dat, 'impossible')=TRUE; return(dat)}
	# # afreq = attr(marker, 'afreq')
	# # impossible = FALSE
	# # sex = x$pedigree[, 'SEX']
	# # dat = lapply(1:x$nInd, function(i) {
		# # h = glist[[i]]
		# # if(i %in% x$founder) {
			# # prob = switch(sex[i], afreq[h], afreq[h[1,]] * afreq[h[2,]] * ((h[1,] != h[2,]) + 1))
			# # if(sum(prob)==0) impossible = TRUE
		# # }
		# # else  prob = rep.int(1, length(h)/sex[i])	
		# # list(hap = h, prob = as.numeric(prob))
	# # })
	# # attr(dat, 'impossible') = impossible
	# # dat
# # }


# # .likel <- function(x, marker1, eliminate=0, startdata=NULL, logbase=NULL) {
	
	# # if(x$hasLoops) stop("Unbroken loops in pedigree.")
	# # #nInd=x$nInd; ped=x$pedigree; chrom=x$model$chrom

	# # if(is.null(startdata)) startdata = .startdata_M(x, marker1, eliminate=eliminate)
	# # dat = startdata
	# # if (attr(dat, 'impossible')) return(ifelse(is.numeric(logbase), -Inf, 0))
	# # uninf =.dropIndivs(x, marker1, peel=FALSE, internal.ids=TRUE)
	
	# # if(is.null(dups <- x$loop_breakers)) {
		# # for (sub in x$subnucs) 	{
			# # pivtype = sub$pivtype
			# # if ((pivtype == 1 && all(unlist(sub[2:3]) %in% uninf)) || (pivtype == 2 && all(unlist(sub[c(1,3)]) %in% uninf))) next #if all offs and the non-piv parent are uninformative
			# # dat = .peel_M(dat, sub, chrom='AUTOSOMAL', SEX=NULL)
			# # if (pivtype > 0 && attr(dat, 'impossible')) return(ifelse(is.numeric(logbase), -Inf, 0))
		# # }
		# # likelihood = dat
	# # }
	# # else {
		# # two2one = function (matr) 1000 * matr[1, ] + matr[2, ]
		# # origs = match(dups[, 1], x$orig.ids); copies = match(dups[, 2], x$orig.ids)

		# # #For each orig, find the indices of its haplos (in orig$hap) that also occur in its copy. Then take cross product of these vectors.
		# # sumover = .my.grid(lapply(seq_along(origs), function(i)	which(two2one(dat[[c(origs[i], 1)]]) %in% two2one(dat[[c(copies[i], 1)]]))), as.list=TRUE)	
		# # likelihood = 0
		# # for(r in sumover) { #r a vector of indices: r[i] gives a column number of the hap matrix of orig[i].
			# # dat1 = dat; attr(dat1, 'impossible') = FALSE
			# # for(i in seq_along(origs)) {
				# # origdat = dat[[origs[i]]]; 
				# # hap = origdat[[1]][, r[i], drop=F]
				# # prob = origdat[[2]][r[i]]; if(sum(prob)==0) print("Loop-loekke: Alle sannsynligheter er null. Magnus lurer paa om dette gir feilmelding.")
				# # dat1[[origs[i]]] = list(hap = hap, prob = prob)
				# # dat1[[copies[i]]] = list(hap = hap, prob = 1)
			# # }
			
			# # for (sub in x$subnucs) 	{
				# # pivtype = sub$pivtype
				# # if ((pivtype == 1 && all(unlist(sub[2:3]) %in% uninf)) || (pivtype == 2 && all(unlist(sub[c(1,3)]) %in% uninf))) next #if all offs and the non-piv parent are uninformative
				# # dat1 = .peel_M(dat1, sub, chrom='AUTOSOMAL', SEX=NULL)
				# # if (pivtype > 0 && attr(dat1, 'impossible')) {break} #if impossible data - break out of ES-algorithm and go to next r in sumover.
				# # if (pivtype == 0) likelihood = likelihood + dat1
			# # }
		# # }
	# # }
	# # if (is.numeric(logbase)) log(likelihood, logbase) else likelihood
# # }



# # .peel_M_X <- function(dat, sub, SEX) { #cat("Parents:", far, mor,"\n")
	# # far = sub[['father']]; mor = sub[['mother']]; offs = sub[['offspring']]; piv = sub[['pivot']]; pivtype = sub[['pivtype']]
	# # if(pivtype==3) offs = offs[offs != piv]
	# # farh = dat[[c(far, 1)]]; morh = dat[[c(mor, 1)]]
	# # likel = dat[[c(far, 2)]] %*% t.default(dat[[c(mor, 2)]])
	# # dims = dim(likel);	fa_len = dims[1L];  mo_len = dims[2L]

	# # for (b in offs) {
		# # datb = dat[[b]]; bh = datb[[1]]; bp = datb[[2]]; bl = length(bp)
		# # switch(SEX[b],
		# # { trans_mats = .trans_M(morh, bh)
		  # # mm = rep(.colSums(trans_mats * bp, bl, mo_len), each=fa_len); 
		# # },
		# # { trans_pats = unlist(lapply(farh, function(fh) as.numeric(fh == bh[1, ])))
		  # # trans_mats = .trans_M(morh, bh[2, ])
		  # # dim(trans_mats) = c(bl, mo_len); trans_mats_rep = as.numeric(do.call(rbind, rep(list(trans_mats), fa_len))) #TODO: improve with 'rep'?
		  # # mm = .colSums((trans_pats * bp) * trans_mats_rep, bl, fa_len*mo_len)
		# # })
		# # likel = likel * mm 
	# # }
	# # if(pivtype==0) return(sum(likel))
		
	# # res = switch(pivtype, .rowSums(likel, fa_len, mo_len), .colSums(likel, fa_len, mo_len),
		# # {	pivh = dat[[c(piv, 1)]]; pivp = dat[[c(piv, 2)]]; pi_len = length(pivp)
			# # switch(SEX[piv],
			# # {	trans_mats = .trans_M(morh, pivh)
				# # T = rep(trans_mats, each=fa_len)
			# # }, 
			# # {	T = numeric(fa_len * mo_len * pi_len);	dim(T) = c(fa_len, mo_len, pi_len)
				# # trans_mats = .trans_M(morh, pivh[2, ]); dim(trans_mats) = c(pi_len, mo_len)
				# # for (i in seq_len(fa_len))
					# # T[i,,] = t.default(as.numeric(farh[i] == pivh[1,]) * trans_mats) #TODO:make faster?
			# # })
			# # arr = as.vector(T) * as.vector(likel); dim(arr) = dim(T)
			# # res = .colSums(arr, fa_len*mo_len, pi_len) #sum for each entry of haps[[piv]]
			# # res * pivp
		# # })
	# # pivhap_update = switch(SEX[piv], dat[[c(piv, 1)]][res>0], dat[[c(piv, 1)]][, res>0, drop=F])
	# # dat[[piv]] = list(hap = pivhap_update, prob = res[res>0])
	# # if(sum(res)==0) attr(dat, 'impossible') = TRUE
	# # return(dat)
# # }


# # .likel_X <- function(x, marker1, eliminate=0, startdata=NULL, logbase=NULL) { #if(exists("ii")) ii<<-ii+1
	
	# # if(x$hasLoops) stop("Unbroken loops in pedigree.")
	# # SEX=x$pedigree[, 'SEX']

	# # if(is.null(startdata)) startdata = .startdata_M_X(x, marker1, eliminate=eliminate)
	# # dat = startdata
	# # if (attr(dat, 'impossible')) return(ifelse(is.numeric(logbase), -Inf, 0))
	# # uninf =.dropIndivs(x, marker1, peel=FALSE, internal.ids=TRUE)
	
	# # if(is.null(lb <- x$loop_breakers)) {
		# # for (sub in x$subnucs) 	{
			# # pivtype = sub$pivtype
			# # if ((pivtype == 1 && all(unlist(sub[2:3]) %in% uninf)) || (pivtype == 2 && all(unlist(sub[c(1,3)]) %in% uninf))) next #if all offs and the non-piv parent are uninformative
			# # dat = .peel_M_X(dat, sub, SEX=SEX)
			# # if (pivtype > 0 && attr(dat, 'impossible')) return(ifelse(is.numeric(logbase), -Inf, 0))
		# # }
		# # likelihood = dat
	# # }
	# # else {
		# # two2one = function (matr) 1000 * matr[1, ] + matr[2, ]
		# # origs = match(lb[, 1], x$orig.ids); copies = match(lb[, 2], x$orig.ids)

		# # #For each orig, find the indices of its haplos (in orig$hap) that also occur in its copy. Then take cross product of these vectors.
		# # sumover = .my.grid(lapply(seq_along(origs), function(i)	{ ori = dat[[c(origs[i], 1)]]; cop = dat[[c(copies[i], 1)]]; if(SEX[origs[i]]==2) {ori=two2one(ori); cop=two2one(cop)}; seq_along(ori)[ori %in% cop]}), as.list=TRUE)	
		# # likelihood = 0
		# # for(r in sumover) { #r a vector of indices: r[i] gives a column number of the hap matrix of orig[i].
			# # dat1 = dat; attr(dat1, 'impossible') = FALSE
			# # for(i in seq_along(origs)) {
				# # orig.int = origs[i]; copy.int = copies[i]
				# # origdat = dat[[orig.int]]; 
				# # hap = switch(SEX[orig.int], origdat[[1]][r[i]], origdat[[1]][, r[i], drop=F])
				# # prob = origdat[[2]][r[i]]; if(sum(prob)==0) print("All probabilities are zero. Will this cause an error?")
				# # dat1[[orig.int]] = list(hap = hap, prob = prob)
				# # dat1[[copy.int]] = list(hap = hap, prob = 1)
			# # }
			# # for (sub in x$subnucs) 	{
				# # pivtype = sub$pivtype
				# # if ((pivtype == 1 && all(unlist(sub[2:3]) %in% uninf)) || (pivtype == 2 && all(unlist(sub[c(1,3)]) %in% uninf))) next #if all offs and the non-piv parent are uninformative
				# # dat1 = .peel_M_X(dat1, sub, SEX=SEX)
				# # if (pivtype > 0 && attr(dat1, 'impossible')) break #if impossible data - break out of ES-algorithm and go to next r in sumover.
				# # if (pivtype == 0) likelihood = likelihood + dat1
			# # }
		# # }
	# # }
	# # if (is.numeric(logbase)) log(likelihood, logbase) else likelihood
# # }



# .startdata_MD_X = function(x, marker) {
	# startprob_X <- function(h, model, afreq, sex, aff, founder) {
 		# if(sex==1) {
			# mat = h #vector
			# d.no = as.numeric(mat<0)
			# prob = switch(aff+1, rep.int(1, length(d.no)), (1 - model$penetrances$male)[d.no+1], model$penetrances$male[d.no+1])
			# if (founder)  
				# prob = prob * afreq[abs(mat)] * c(1-model$dfreq, model$dfreq)[d.no+1]
		# } else {
			# pat = h[1,]; mat = h[2,]
			# al1 = abs(pat); al2 = abs(mat)
			# d.no = (pat<0) + (mat<0)
			# prob = switch(aff+1, rep.int(1, length(d.no)), (1 - model$penetrances$female)[d.no+1], model$penetrances$female[d.no+1])
			# if (founder)
				# prob = prob * afreq[al1] * afreq[al2] * ((al1 != al2) + 1) * model$dfreq^d.no * (1-model$dfreq)^(2-d.no)
		# }
		# as.numeric(prob)
  	# }
	
	# ## MAIN ###
	# marker = .reduce_alleles(marker)
	# glist = .build_genolist_X(x, marker)
	# if (attr(glist, 'impossible')) {
		# dat=list(); attr(dat, 'impossible')=TRUE; return(dat)
	# }

	# dfact_fem = cbind(c(-1,-1),c(-1,1),c(1,-1),c(1,1)) #D=-1; N=1
	# SEX = x$pedigree[,'SEX']; AFF = x$pedigree[,'AFF']; FOU = (1:x$nInd) %in% x$founder
	# impossible = FALSE
	# dat = lapply(1:x$nInd, function(i) {
		# switch(SEX[i], 
		# { hap = c(glist[[i]], -glist[[i]])
		  # prob = startprob_X(hap, model=x$model, afreq=attr(marker, 'afreq'), sex=1, aff=AFF[i], founder=FOU[i])
		  # list(hap = hap[prob > 0], prob = prob[prob > 0])},
		# { gl = ncol(glist[[i]])  
		  # hap = glist[[i]][, rep(1:gl, each=4), drop=F ] *  dfact_fem[, rep(1:4, times=gl), drop=F]
		  # prob = startprob_X(hap, model=x$model, afreq=attr(marker, 'afreq'), sex=2, aff=AFF[i], founder=FOU[i])
		  # if (all(prob==0)) impossible = TRUE
		  # list(hap = hap[, prob > 0, drop=F], prob = as.numeric(prob[prob > 0]))
		# })
	# })
	# attr(dat, 'impossible') = impossible
	# dat
# }


	
# .likelihood_multi_X <- function(x, theta, dat, logbase=NULL) {
	
	# if(x$hasLoops) stop("Unbroken loops in pedigree.")
	# nInd=x[['nInd']]; ped=x[['pedigree']]; SEX=ped[, 'SEX']; chrom=x[['model']][['chrom']]
	
	# if(is.null(dups <- x$loop_breakers)) {
		# for (sub in x$subnucs) 	{
			# dat = .peel_MD_X(dat, sub, theta, SEX)
			# if (sub$pivtype > 0 && attr(dat, 'impossible')) return(ifelse(is.numeric(logbase), -Inf, 0))
		# }
		# likelihood = dat
	# }
	# else {
		# if (attr(dat, 'impossible')) return(ifelse(is.numeric(logbase), -Inf, 0))
		# two2one = function (matr) 1000 * matr[1, ] + matr[2, ]
		# origs = match(dups[, 1], x$orig.ids); copies = match(dups[, 2], x$orig.ids)

		# #For each orig, find the indices of its haplos (in orig$hap) that also occur in its copy. Then take cross product of these vectors.
		# sumover = .my.grid(lapply(seq_along(origs), function(i)
			# switch(SEX[origs[i]], 
			 # which(dat[[c(origs[i], 1)]] %in% dat[[c(copies[i], 1)]]),
			 # which(two2one(dat[[c(origs[i], 1)]]) %in% two2one(dat[[c(copies[i], 1)]]))
			# )), as.list=TRUE)
		
		# likelihood = 0
		# for(r in sumover) { #r a vector of indices: r[i] gives a column number of the hap matrix of orig[i].
			# dat1 = dat; 
			# for(i in seq_along(origs)) {
				# origdat = dat[[origs[i]]]; 
				# hap = switch(SEX[origs[i]], origdat[[1]][r[i]], origdat[[1]][, r[i], drop=F])
				# dat1[[origs[i]]] = list(hap = hap, prob = origdat[[2]][r[i]])
				# dat1[[copies[i]]] = list(hap = hap, prob = 1)
			# }
			# for (sub in x$subnucs) 	{
				# dat1 = .peel_MD_X(dat1, sub, theta, SEX)
				# if (sub$pivtype > 0 && attr(dat1, 'impossible')) break #if impossible data - break out of ES-algorithm and go to next r in sumover.
				# if (is.numeric(dat1)) {likelihood = likelihood + dat1; break}
			# }
		# }
	# }
	# if (is.numeric(logbase)) log(likelihood, logbase) else likelihood
# }


# .peel_MD_X <- function(dat, sub, theta, SEX) { 
	# far = sub[['father']]; mor = sub[['mother']]; offs = sub[['offspring']]; piv = sub[['pivot']]; pivtype = sub[['pivtype']]
	# if(pivtype==3) offs = offs[offs != piv]
	
	# farh = dat[[c(far, 1)]]; morh = dat[[c(mor, 1)]]
	# likel = dat[[c(far, 2)]] %*% t.default(dat[[c(mor, 2)]])
	# dims = dim(likel);	fa_len = dims[1L];  mo_len = dims[2L];
	# mm_init = numeric(length(likel)); dim(mm_init) = dim(likel)

	# for (b in offs) {
		# mm = mm_init
		# switch(SEX[b],
			# {for (j in seq_len(mo_len)) 
				# mm[, j] = .trans_MD(morh[, j], dat[[c(b,1)]], theta) %*% dat[[c(b,2)]]
			# }, 
			# {for (j in seq_len(mo_len)) {
				# transmother = .trans_MD(morh[, j], dat[[c(b,1)]][2,], theta)
				# for (i in seq_len(fa_len))
					# mm[i,j] = (as.numeric(farh[i] == dat[[c(b,1)]][1,]) * transmother) %*% dat[[c(b,2)]]	
			# }}
		# )
		# likel = likel * mm
	# }
	# if(pivtype==0) return(sum(likel))
	
	# res = switch(pivtype, .rowSums(likel, fa_len, mo_len), .colSums(likel, fa_len, mo_len),
	# {
		# pivh = dat[[c(piv, 1)]]; pivp = dat[[c(piv, 2)]]; pi_len = length(pivp)
		# T = numeric(fa_len * mo_len * pi_len);	dim(T) = c(fa_len, mo_len, pi_len)
		# switch(SEX[piv],
		# {for (j in seq_len(mo_len))
			# transmother = .trans_MD(morh[, j], pivh, theta)
			# for(i in seq_len(fa_len)) T[i,j,] = transmother
		# }, 
		# {for (j in seq_len(mo_len)) {
			# transmother = .trans_MD(morh[, j], pivh[2,], theta)
			# for (i in seq_len(fa_len))
				# T[i,j,] = (as.numeric(farh[i] == pivh[1,]) * transmother)
		# }})

		# arr = as.vector(T) * as.vector(likel); dim(arr) = dim(T)
		# res = .colSums(arr, fa_len*mo_len, pi_len) #sum for each entry of haps[[piv]]
		# res * pivp
	# })
	# pivhap_update = switch(SEX[piv], dat[[c(piv, 1)]][res>0], dat[[c(piv, 1)]][, res>0, drop=F])
	# dat[[piv]] = list(hap = pivhap_update, prob = res[res>0])
	# if(all(res==0)) attr(dat, 'impossible') = TRUE
	# return(dat)
# }


# .likelihood_multi <- function(x, theta, dat, logbase=NULL) {
	
	# if(x$hasLoops) stop("Unbroken loops in pedigree.")
	# nInd=x[['nInd']]; ped=x[['pedigree']]; chrom=x[['model']][['chrom']]
	# if (attr(dat, 'impossible')) return(ifelse(is.numeric(logbase), -Inf, 0))
	
	# if(is.null(dups <- x$loop_breakers)) {
		# for (sub in x$subnucs) 	{
			# dat = .peel_MD(dat, sub, theta, chrom='AUTOSOMAL', SEX=NULL)
			# if (sub$pivtype > 0 && attr(dat, 'impossible')) return(ifelse(is.numeric(logbase), -Inf, 0))
		# }
		# likelihood = dat
	# }
	# else {
		# two2one = function (matr) 1000 * matr[1, ] + matr[2, ]
		# origs = match(dups[, 1], x$orig.ids); copies = match(dups[, 2], x$orig.ids)

		# #For each orig, find the indices of its haplos (in orig$hap) that also occur in its copy. Then take cross product of these vectors.
		# sumover = .my.grid(lapply(seq_along(origs), function(i)	which(two2one(dat[[c(origs[i], 1)]]) %in% two2one(dat[[c(copies[i], 1)]]))), as.list=TRUE)
		# likelihood = 0
		# for(r in sumover) { #r a vector of indices: r[i] gives a column number of the hap matrix of orig[i].
			# dat1 = dat; 
			# for(i in seq_along(origs)) {
				# origdat = dat[[origs[i]]]; 
				# hap = origdat[[1]][, r[i], drop=F]
				# dat1[[origs[i]]] = list(hap = hap, prob = origdat[[2]][r[i]])
				# dat1[[copies[i]]] = list(hap = hap, prob = 1)
			# }
			# for (sub in x$subnucs) 	{
				# dat1 = .peel_MD(dat1, sub, theta, chrom='AUTOSOMAL', SEX=NULL)
				# if (sub$pivtype > 0 && attr(dat1, 'impossible')) break #if impossible data - break out of ES-algorithm and go to next r in sumover.
				# if (is.numeric(dat1)) {likelihood = likelihood + dat1; break}
			# }
		# }
	# }
	# if (is.numeric(logbase)) log(likelihood, logbase) else likelihood
# }


# .eliminate_OLD_X = function(x, genolist) {
	# SEX = x$pedigree[, 'SEX']
	# for (sub in x$subnucs) {
		# mor = genolist[[ sub[[3]] ]]; barnlist=sub[-c(1:3)]
		# far1 = genolist[[ sub[[2]] ]] #already a vector, se genolist.
		# fl = length(far1)
		# mor_ext1 = rep(mor[1,], each=fl); mor_ext2 = rep(mor[2,], each=fl); ml = ncol(mor);
# #	
		# comp <- lapply(barnlist, function(barn) {
			# b = genolist[[barn]]
			# compv = switch(SEX[barn],
			# {  b2 = rep(b, each=fl*ml); (mor_ext1 == b2 | mor_ext2 == b2) },
			# {  b1 = rep(b[1,], each=fl*ml); b2 = rep(b[2,], each=fl*ml)
			  # (far1 == b1) & (mor_ext1 == b2 | mor_ext2 == b2) 
			# })
			# dim(compv) = c(fl, ml, length(b)/SEX[barn])
			# compv
		# })			
		# parent_keep = rep.int(T, fl*ml)
		# for(arr in comp) 
			# parent_keep = parent_keep * rowSums(arr, dims=2) #sums over dimension 3, i.e. parent_keep goes to 0 if there is an offsp with no compat genotypes. Faster than apply(,3,any)
		# far_keep = rowSums(parent_keep)>0
		# mor_keep = colSums(parent_keep)>0
		# genolist[[ sub[[2]] ]] <- far1[far_keep]
		# genolist[[ sub[[3]] ]] <- mor[, mor_keep, drop=F]
		# offs_keep <- lapply(comp, function(arr) colSums(arr[far_keep, mor_keep,,drop=F], dims=2)>0) #sums over dimensions 1:2
		# for (bno in seq_along(barnlist)) {
			# barn = barnlist[bno]
			# barngeno = genolist[[barn]]
			# keep = offs_keep[[bno]]
			# genolist[[barn]] <- switch(SEX[barn], barngeno[keep], barngeno[ , keep, drop=F])
		# }
	# }
	# if(any(unlist(lapply(genolist[sub[-1]], length))) == 0) {attr(genolist, 'impossible') = TRUE; return(genolist)}
	# genolist
# }


# .eliminate_OLD = function(x, genolist) {
	# for (sub in x$subnucs) {
		# far = genolist[[ sub[['father']] ]]; mor = genolist[[ sub[['mother']] ]]; barnlist=sub[['offspring']]
		# far1 = far[1,]; far2 = far[2,]; 
		# fl = length(far1); ml = ncol(mor)
		# mor_ext1 = rep(mor[1,], each=fl); mor_ext2 = rep(mor[2,], each=fl)
		
		# comp = lapply( genolist[barnlist], function(b) {
			# b1 = rep(b[1,], each=fl*ml); b2 = rep(b[2,], each=fl*ml)
			# compv = (far1 == b1 | far2 == b1) & (mor_ext1 == b2 | mor_ext2 == b2) #far and mor_ext recycle
			# dim(compv) = c(fl, ml, ncol(b))
			# compv
		# })			
		# parent_keep = rep.int(T, fl*ml)
		# for(arr in comp) 
			# parent_keep = parent_keep * .rowSums(arr, fl*ml, dim(arr)[3]) #rowSums(arr, dims=2) #sums over dimension 3, i.e. parent_keep goes to 0 if there is an offsp with no compat genotypes. Faster than apply(,3,any)
		# far_keep = .rowSums(parent_keep, fl, ml)>0
		# mor_keep = .colSums(parent_keep, fl, ml)>0
		# genolist[[ sub[['father']] ]] <- far[, far_keep, drop=F]
		# genolist[[ sub[['mother']] ]] <- mor[, mor_keep, drop=F]
		# offs_keep <- lapply(comp, function(arr) {arr_red = arr[far_keep, mor_keep,,drop=F]; dn = dim(arr_red); .colSums(arr_red, dn[1L]*dn[2L], dn[3L])>0}) #sums over dimensions 1:2
		# for (bno in seq_along(barnlist)) genolist[[barnlist[bno]]] <- genolist[[barnlist[bno]]][ , offs_keep[[bno]], drop=F]
		# if(any(unlist(lapply(genolist[unlist(sub[1:3])], length))) == 0) {attr(genolist, 'impossible') = TRUE; return(genolist)}
	# }
	# genolist
# }



# .likelihood_MM <- function(x, theta, marker1, marker2, startdata=NULL, eliminate=0, logbase=NULL) {
	
	# if(x$hasLoops) stop("Unbroken loops in pedigree.")
	# #nInd=x[['nInd']]; ped=x[['pedigree']]; chrom=x[['model']][['chrom']]
	
	# if(is.null(startdata)) startdata = .startdata_MM(x, marker1, marker2, eliminate=eliminate)
	# dat = startdata
	
	# if(is.null(dups <- x$loop_breakers)) {
		# for (sub in x$subnucs) 	{
			# if (any(unlist(lapply(dat, function(dati) length(dati$prob) == 0)))) return(ifelse(is.numeric(logbase), -Inf, 0))
			# dat = .peel_MM(dat, sub, theta)
		# }
		# likelihood = dat
	# }
	# else {
		# two2one = function (matr) 1000 * matr[1, ] + matr[2, ]
		# origs = match(dups[, 1], x$orig.ids); copies = match(dups[, 2], x$orig.ids)

		# #For each orig, find the indices of its haplos (in orig$hap) that also occur in its copy. Then take cross product of these vectors.
		# sumover = .my.grid(lapply(seq_along(origs), function(i) which(two2one(dat[[c(origs[i], 1)]]) %in% two2one(dat[[c(copies[i], 1)]]))), as.list=TRUE)
		# likelihood = 0
		# for(r in sumover) { #r a vector of indices: r[i] gives a column number of the hap matrix of orig[i].
			# dat1 = dat; 
			# for(i in seq_along(origs)) {
				# origdat = dat[[origs[i]]]; 
				# hap = origdat[[1]][, r[i], drop=F]
				# dat1[[origs[i]]] = list(hap = hap, prob = origdat[[2]][r[i]])
				# dat1[[copies[i]]] = list(hap = hap, prob = 1)
			# }
			# for (sub in x$subnucs) 	{
				# dat1 = .peel_MM(dat1, sub, theta)
				# if (is.numeric(dat1)) {likelihood = likelihood + dat1; break}
				# if (any(unlist(lapply(dat1, function(dat1i) length(dat1i$prob) == 0)))) break #if impossible data - break out of ES-algorithm and go to next r in sumover.
			# }
		# }
	# }
	# if (is.numeric(logbase)) log(likelihood, logbase) else likelihood
# }

exclusionPower = function(ped_claim, ped_true, ids, alleles, afreq=NULL, known_genotypes=list(), loop_breakers=NULL, Xchrom=FALSE) {
	if(class(ped_claim)=='linkdat') 
      ped_claim = list(ped_claim)
   if(class(ped_true)=='linkdat') 
      ped_true = list(ped_true)
   
   singletons_claim = sapply(ped_claim, function(x) is.numeric(x) && length(x)==1)
   singletons_true = sapply(ped_true, function(x) is.numeric(x) && length(x)==1)  
   if(any(singletons_claim))
      ped_claim[singletons_claim] = lapply(ped_claim[singletons_claim], function(id) {
         lables = id*100000+c(1,2,3) #one of the parents will become 'id'
         for(ped in ped_true[!singletons_true]) #check the other pedigrees
            if(id %in% ped$orig.ids) {
               sex = paramlink:::.getSex(ped, id)
               lables[sex] = id
               break
            }
         relabel(nuclearPed(1), lables)})
   if(!all(sapply(ped_claim, class) == 'linkdat')) stop("Something is wrong with the 'ped_claim' input") 
   if(any(singletons_true))
      ped_true[singletons_true] = lapply(ped_true[singletons_true], function(id) {
         lables = id*100000+c(1,2,3) #one of the parents will become 'id'
         for(ped in ped_claim[!singletons_claim]) #check the other pedigrees 
            if(id %in% ped$orig.ids) {
               sex = paramlink:::.getSex(ped, id)
               lables[sex] = id
               break
            }
         relabel(nuclearPed(1), lables)})
   if(!all(sapply(ped_true, class) == 'linkdat')) stop("Something is wrong with the 'ped_true' input") 
   
   ids_claim = lapply(ped_claim, function(x) ids[ids %in% x$orig.ids])
   ids_true = lapply(ped_true, function(x) ids[ids %in% x$orig.ids])
   loops_claim = lapply(ped_claim, function(x) {lb = x$orig.ids[x$orig.ids %in% loop_breakers]; if(length(lb)==0) lb = NULL; lb})
   loops_true = lapply(ped_true, function(x) {lb = x$orig.ids[x$orig.ids %in% loop_breakers]; if(length(lb)==0) lb = NULL; lb})
   
   N_claim = length(ped_claim)
   N_true = length(ped_true)
   
   if(length(alleles)==1) alleles=1:alleles
   if(Xchrom) chrom = 23 else chrom = NA
   partial_claim = lapply(1:N_claim, function(i) {
      x = ped_claim[[i]]
      m = marker(x, alleles=alleles, afreq=afreq, chrom=chrom)
      for(tup in known_genotypes) 
         if(tup[1] %in% x$orig.ids)
            m = modifyMarker(x, m, ids=tup[1], genotype=tup[2:3])
      m})
   partial_true = lapply(1:N_true, function(i) {
      x = ped_true[[i]]
      m = marker(x, alleles=alleles, afreq=afreq, chrom=chrom)
      for(tup in known_genotypes) 
         if(tup[1] %in% x$orig.ids)
            m = modifyMarker(x, m, ids=tup[1], genotype=tup[2:3])
      m})
   has_genotypes = length(known_genotypes)>0
   oldmfrow = par()[['mfrow']]
   par(mfrow=c(1, N_claim + N_true))
   for(i in 1:N_claim) {
      x = ped_claim[[i]]
      mm = if(has_genotypes) partial_claim[[i]] else NULL
      tit = if(N_claim==1) 'Claim' else paste("Claim, part", i)
      cols = ifelse(x$orig.ids %in% ids_claim[[i]], 2, 1)
      lables = sapply(x$orig.ids, function(i) if(i>100000) 'dummy' else i)
      plot(x, marker=mm, title=tit, col=cols, id.labels=lables)
   }
   for(i in 1:N_true) {
      x = ped_true[[i]]
      mm = if(has_genotypes) partial_true[[i]] else NULL
      tit = if(N_true==1) 'True' else paste("True, part", i)
      cols = ifelse(x$orig.ids %in% ids_true[[i]], 2, 1)
      lables = sapply(x$orig.ids, function(i) if(i>100000) 'dummy' else i)
      plot(x, marker=mm, title=tit, col=cols, id.labels=lables)
   }
   par(mfrow=oldmfrow)
   p.g = Reduce('%o%', lapply(1:N_true, function(i) 
      oneMarkerDistribution(ped_true[[i]], ids=ids_true[[i]], partialmarker=partial_true[[i]], loop_breakers=loops_true[[i]], verbose=F)))
   I.g = Reduce('%o%', lapply(1:N_claim, function(i) 
      oneMarkerDistribution(ped_claim[[i]], ids=ids_claim[[i]], partialmarker=partial_claim[[i]], loop_breakers=loops_claim[[i]], verbose=F)==0))
   sum(p.g*I.g)
}

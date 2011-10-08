nuclearPed <- function(boys, girls=0) {
	stopifnot(boys+girls>0)
	p=cbind(ID=1:(2+boys+girls), FID=c(0, 0, rep(1, girls+boys)), MID=c(0, 0, rep(2, girls+boys)),
			SEX=c(1, 2, rep(1, boys), rep(2, girls)), AFF=1)
	linkdat(p, verbose=FALSE)
}

cousinPed <- function(degree) {
	stopifnot(degree>=0)
	if(degree==0) return(nuclearPed(boys=1, girls=1))
	p = cbind(ID=1:4, FID=c(0,0,1,1), MID=c(0,0,2,2), SEX=c(1,2,1,1), AFF=1)
	for (n in 1:degree)
		p = rbind(p, c(4*n+1,0,0,2,1), c(4*n+2,0,0,2,1), c(4*n+3,4*n-1,4*n+1,1,1), c(4*n+4,4*n,4*n+2,1,1))
	p[nrow(p), 'SEX'] = 2
	linkdat(p, verbose=FALSE)
}

halfCousinPed <- function(degree) {
	stopifnot(degree>=0)
	if(degree==0) p = cbind(ID=1:5, FID=c(0,0,0,1,1), MID=c(0,0,0,2,3), SEX=c(1,2,2,1,2), AFF=1)
	else {
		p = cbind(ID=1:3, FID=c(0,0,0), MID=c(0,0,0), SEX=c(1,2,2), AFF=1)
		for (n in seq_len(degree))
			p = rbind(p, c(4*n,4*n-4,4*n-3,1,1), c(4*n+1,0,0,2,1), c(4*n+2,4*n-2,4*n-1,1,1), c(4*n+3,0,0,2,1)) #add 1 generation: son in line 1, his wife, son in line 2, his wife
		dd=degree+1
		p = rbind(p, c(4*dd,4*dd-4,4*dd-3,1,1), c(4*dd+1,4*dd-2,4*dd-1,2,1)) #last generation - one boy, one girl.
		p[4,c(2,3)] = c(1,2); p[6,c(2,3)] = c(1,3) 
	}
	linkdat(p, verbose=FALSE)
}

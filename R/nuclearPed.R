nuclearPed <- function(boys, girls=0) {
	stopifnot(boys+girls>0)
	p=cbind(ID=1:(2+boys+girls), FID=c(0, 0, rep(1, girls+boys)), MID=c(0, 0, rep(2, girls+boys)),
			SEX=c(1, 2, rep(1, boys), rep(2, girls)), AFF=1)
	linkdat(p, verbose=FALSE)
}

cousinPed <- function(degree=1) {
	stopifnot(degree>0)
	p = cbind(ID=1:4, FID=c(0,0,1,1), MID=c(0,0,2,2), SEX=c(1,2,1,1), AFF=1)
	for (n in 1:degree)
		p = rbind(p, c(4*n+1,0,0,2,1), c(4*n+2,0,0,2,1), c(4*n+3,4*n-1,4*n+1,1,1), c(4*n+4,4*n,4*n+2,1,1))
	linkdat(p, verbose=FALSE)
}

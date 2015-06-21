
dbRDA.D <- function(D, X, nperm=999, option=3, compute.eig=FALSE, coord=FALSE, rda.coord=2, positive.RDA.values=FALSE)
{
	D <- as.matrix(D)
	X <- as.matrix(X)
	n <- nrow(D)
	epsilon <- sqrt(.Machine$double.eps)
#
# Gower centring, matrix formula. Legendre & Legendre (2012), equation 9.42
	One <- matrix(1,n,n)
	mat <- diag(n) - One/n
	G <- -0.5 * mat %*% (D^2) %*% mat
	SSY <- sum(diag(G))
	# LCBD <- diag(G)
#
# Principal coordinate analysis after eigenvalue decomposition of D
	if(compute.eig) {
		eig <- eigen(G, symmetric=TRUE)
		values <- eig$values     # All eigenvalues
		vectors <- eig$vectors   # All eigenvectors, scaled to lengths 1
		if(coord) {
			select <- which(values > epsilon)
			princ.coord <- vectors[,select] %*% diag(sqrt(values[select]))
			} else { princ.coord <- NA }
		} else {
		values <- princ.coord <- NA
		}
#
# Compute projector matrix H ("hat" matrix in the statistical literature)
	X.c <- scale(X, center=TRUE, scale=FALSE)   # Centre matrix X
	m <- qr(X.c, tol=1e-6)$rank                 # m = rank of X.c
	cat("Rank of X centred =",m,"\n")
	if(m==1) { 
		H <- (X.c[,1] %*% t(X.c[,1]))/((t(X.c[,1]) %*% X.c[,1])[1,1]) 
		} else {
		if(option<3) {  
			# if(det(t(X.c)%*%X.c)<epsilon) stop ('Collinearity detected in X')
			if(m < ncol(X.c)) stop ('Collinearity detected in X')
			H <- X.c %*% solve(t(X.c) %*% X.c) %*% t(X.c)
			#
			# option=3: compute projector H from orthogonalized X; no inversion
			} else {
			X.eig <- eigen(cov(X.c))
			k <- length(which(X.eig$values > epsilon))
                        X.ortho <- X.c %*% X.eig$vectors[,1:k]  # F matrix of PCA
			XprX <- t(X.ortho) %*% X.ortho
			H <- X.ortho %*% diag(diag(XprX)^(-1)) %*% t(X.ortho)
			}
	}
#
# Compute the F statistic: McArdle & Anderson (2001), equation 4 modified
	HGH <- H %*% G %*% H
	SSYhat <- sum(diag(HGH))
	#
	if(option==1) {
		I.minus.H <- diag(n) - H
		den1 <- sum(diag(I.minus.H %*% G %*% I.minus.H))
		F <- SSYhat/den1     # F statistic without the degrees of freedom
		Rsquare <- F/(F+1)
	} else {
		F <- SSYhat/(SSY-SSYhat)  # F statistic without the degrees of freedom
		Rsquare <- SSYhat/SSY     # or equivalent: Rsquare <- F/(F+1)
	}
	RsqAdj <- 1-((1-Rsquare)*(n-1)/(n-1-m))
#
# Permutation test of F
	if(nperm > 0) {
		nGE=1
		for(i in 1:nperm) {
			order <- sample(n)
			Gperm <- G[order, order]
			H.Gperm.H <- H %*% Gperm %*% H
			SSYhat.perm <- sum(diag(H.Gperm.H))
			#
			if(option==1) {
				den <- sum(diag(I.minus.H %*% Gperm %*% I.minus.H))
				F.perm <- SSYhat.perm/den
			} else {
				F.perm <- SSYhat.perm/(SSY-SSYhat.perm)
			}
			if(F.perm >= F) nGE=nGE+1
			}
		P.perm <- nGE/(nperm+1)
		} else { P.perm <- NA }
#
# Compute RDA ordination coordinates
	if(rda.coord > 0) {
		HGH.eig <- eigen(HGH, symmetric=TRUE)
		# kk <- length(which(HGH.eig$values > epsilon))
		RDA.values <- HGH.eig$values
		rel.eig <- RDA.values/SSY
		cum.eig <- cumsum(rel.eig) 
		kk <- length(which(rel.eig > epsilon))
		if(positive.RDA.values) {
			RDA.values <- RDA.values[1:kk]
			rel.eig <- rel.eig[1:kk]
			cum.eig <- cum.eig[1:kk]
			}
		k <- min(rda.coord, kk)
		if(k >= 2) {
		RDA.coord <-sweep(HGH.eig$vectors[,1:k],2,sqrt(RDA.values[1:k]),FUN="*")
			} else {
			RDA.coord <- NA
			cat("k =",k," -- Fewer than two RDA eigenvalues > 0\n")
			}
		} else { RDA.values <- rel.eig <- cum.eig <- RDA.coord <- NA }
#
list(F=F*(n-m-1)/m, Rsquare=c(Rsquare,RsqAdj), P.perm=P.perm, SS.total=SSY, PCoA.values=values, PCoA.vectors=princ.coord, RDA.values=RDA.values/(n-1), RDA.rel.values=rel.eig, RDA.cum.values=cum.eig, RDA.coord=RDA.coord)
}

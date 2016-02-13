RV <-
function(X, Y) {
	#CC# General checks
	if (nrow(X) != nrow(Y)) stop("'X' needs to have the same number of rows as 'Y'")
	if (nrow(X) == 1) stop("Impossible to calculate RV using 1 object")
	
	Y <- scale(Y, scale = FALSE)
	X <- scale(X, scale = FALSE)
	
	XXt <- tcrossprod(X)
	YYt <- tcrossprod(Y)
	
	rv <- sum(diag(XXt %*% YYt))/(sum(diag(XXt %*% XXt)) * sum(diag(YYt %*% YYt)))^0.5
	
	return(rv)
}

"coef.cca" <-
function (object, ...) 
{
	Q <- object$CCA$QR
	u <- object$CCA$u
	u <- sweep(u, 1, sqrt(object$rowsum), "*")
	qr.coef(Q, u)
}


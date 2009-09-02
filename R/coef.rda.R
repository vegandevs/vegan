"coef.rda" <-
function (object, ...) 
{
	Q <- object$CCA$QR
	qr.coef(Q, object$CCA$u)
}


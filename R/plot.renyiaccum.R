`plot.renyiaccum` <-
function (x, what=c("mean", "Qnt 0.025", "Qnt 0.975"), type = "l", ...) 
{
	if (any(what %in% colnames(x[,1,])))
	    x <- x[,,what]
	dm <- dim(x)
	lin <- rep(colnames(x[,1,]), each=dm[1]*dm[2])
	Sites <- rep(1:dm[1], len=prod(dm))
	alp <- factor(rownames(x[1,,]), levels=rownames(x[1,,]))
	alpha <- rep(rep(alp, each=dm[1]), len=prod(dm))
	Diversity <- as.vector(x)
	xyplot(Diversity ~ Sites | alpha, groups=lin, type=type, ...)
}


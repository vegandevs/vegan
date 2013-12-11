`plot.renyiaccum` <-
function (x, what=c("Collector", "mean", "Qnt 0.025", "Qnt 0.975"),
          type = "l", ...) 
{
        what <- what[what %in% dimnames(x)[[3]]]
	if (any(what %in% dimnames(x)[[3]]))
	    x <- x[,,what, drop = FALSE]
	dm <- dim(x)
        dnam <- dimnames(x)
	lin <- rep(dnam[[3]], each=dm[1]*dm[2])
	Sites <- rep(1:dm[1], len=prod(dm))
	alp <- factor(dnam[[2]], levels=dnam[[2]])
	alpha <- rep(rep(alp, each=dm[1]), len=prod(dm))
	Diversity <- as.vector(x)
	xyplot(Diversity ~ Sites | alpha, groups=lin, type=type, ...)
}


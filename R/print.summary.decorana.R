`print.summary.decorana` <-
function (x, head=NA, tail=head, ...) 
{
    digits <- x$digits
    hcat <- function(x, head=head, tail=tail, ...) {
    	   if(!is.na(head) && !is.na(tail) && head + tail + 4 < nrow(x))
    	        x <- rbind(head(x, n=head), "...." = NA, tail(x, n=tail))
    	   printCoefmat(x,  na.print="", ...)
    	}
    if (!is.null(x$spec.scores)) {
        cat("Species scores:\n\n")
        TABLE <- cbind(x$spec.scores, Weights = x$spec.priorweights, 
            Totals = x$spec.totals)
        hcat(TABLE, head=head, tail=tail, digits = digits,  ...)
        cat("\n")
    }
    if (!is.null(x$site.scores)) {
        cat("Site scores:\n\n")
        TABLE <- cbind(x$site.scores, Totals = x$site.totals)
        hcat(TABLE, head=head, tail=tail, digits = digits,  ...)
        cat("\n")
    }
    invisible(x)
}


summary.adipart <-
function(object, digits=3, ...)
{
## internal
p2a <- function(x, y, digits){
    if (length(x) != length(y)) stop("legths differ")
    x <- round(x, digits=digits)
    out <- y2 <- rep(NA, length(x))
    y2[y >= 0.1] <- "ns"
    y2[y < 0.1] <- "."
    y2[y < 0.05] <- "*"
    y2[y < 0.01] <- "**"
    y2[y < 0.001] <- "***"
    for (i in 1:length(x)) out[i] <- paste(as.character(x[i]),y2[i],sep=" ")
    out[is.na(x)] <- ""
    names(out) <- names(x)
    return(out)
}
    x<- object
    test <- attr(x, "times") !=0
    obs.b <- out.b <- round(x$obs$beta[nrow(x$obs$beta):1,], digits)
    obs.a <- out.a <- round(x$obs$alpha[(nrow(x$obs$alpha)-1):1,], digits)
    n <- ncol(out.a)
    Gamma <- round(x$obs$alpha[nrow(x$obs$alpha),],digits)
    if (test) {
        pv.b <- x$exp$beta$p.value[nrow(x$exp$beta$p.value):1,]
        pv.a <- x$exp$alpha$p.value[nrow(x$exp$alpha$p.value):1,]
        if (ncol(x$obs$alpha) != 1) {
            for (i in 1:n) {
                out.b[,i] <- p2a(obs.b[,i], pv.b[,i], digits)
                out.a[,i] <- p2a(obs.a[,i], pv.a[,i], digits)}
        } else {
            out.b <- p2a(obs.b, pv.b, digits)
            out.a <- p2a(obs.a, pv.a, digits)}
        }
        out <- t(data.frame(Gamma, t(out.b), t(out.a)))
        out[is.na(out)] <- ""
        if (ncol(x$obs$alpha) == 1)
            colnames(out) <- colnames(x$obs$alpha)
    output <- list(call=x$input$call, divres=out, times=attr(x, "times"))
    class(output) <- "summary.adipart"
return(output)
}

plot.adipart <-
function(x, rel.yax=NULL, ymax=NULL, p.legend="bottomright", ...)
{
#    x <- object
    if (is.null(rel.yax)) {
        if (attr(x, "design") == "oneway") rel <- TRUE
        if (attr(x, "design") == "twoway") rel <- FALSE
        } else rel <- rel.yax
    if (attr(x, "method") == "trad") {
#        main="Additive diversity partitions"
        if (attr(x, "design") == "oneway")
            xlab <- "Traditional diversity indices"
            else xlab <- "Habitat classes"
        lty <- 1}
    if (attr(x, "method") == "tsallis") {
#        main="Additive Tsallis diversity partitions"
        if (attr(x, "design") == "oneway")
            xlab <- "Scale parameter"
            else xlab <- "Habitat classes"
        lty <- 2}
    if (attr(x, "design") == "oneway") {
        ylab <- if (!rel) "Diversity components" else "Relative diversity components"
        } else {
            if (attr(x, "method") == "tsallis") qspacer <- "q = "
            if (attr(x, "method") == "trad") qspacer <- ""
            spacer <- if (!rel) "D" else "Relative d"
            ylab <- paste(spacer, "iversity partitions (", qspacer, attr(x, "index"), ")", sep="")
            }
    if (!is.null(x$exp)) {
        dat <- list(obs=x$obs$alpha,
            perm=x$exp$alpha$mean,
            pvala=x$exp$alpha$p.value,
            pvalb=x$exp$beta$p.value)
        } else {
        dat <- list(obs=x$obs$alpha,
            perm=x$obs$alpha,
            pvala=x$obs$alpha,
            pvalb=x$obs$alpha)
        dat$pvala[dat$pvala != 1] <- 1
        dat$pvalb[dat$pvalb != 1] <- 1}

    dat$perm <- t(data.frame(t(dat$perm), dat$obs[nrow(dat$obs),]))
    dat$pvala <- t(data.frame(t(dat$pvala), rep(1, ncol(dat$obs))))
    dat$pvalb <- t(data.frame(t(dat$pvalb), rep(1, ncol(dat$obs))))

    if (attr(x, "design") == "twoway") {
        dims <- dim(x$obs$beta)

        dat$obs <- dat$obs[,c(ncol(dat$obs),1:(ncol(dat$obs)-1))]
        dat$perm <- dat$perm[,c(ncol(dat$perm),1:(ncol(dat$perm)-1))]
        dat$pvala <- dat$pvala[,c(ncol(dat$pvala),1:(ncol(dat$pvala)-1))]
        dat$pvalb <- dat$pvalb[,c(ncol(dat$pvalb),1:(ncol(dat$pvalb)-1))]
        dat$obs[dims[1], 2:dims[2]] <- dat$obs[(dims[1] + 1), 2:dims[2]]
        dat$perm[dims[1], 2:dims[2]] <- dat$perm[(dims[1] + 1), 2:dims[2]]
        dat$pvala[dims[1], 2:dims[2]] <- 1
        dat$pvalb[dims[1], 2:dims[2]] <- 1}

    m <- t(dat$obs)
    mp <- t(dat$perm)
    pdot <- t(dat$pvala < 0.05)
    plin <- t(dat$pvalb < 0.05)
    n <- ncol(m)
    lwd <- 5
    col <- "grey"
    col2 <- "black"
    bg <- "white"
    pch <- 21
    nnn <- rownames(m)
    if (length(nnn) == 1) nnn <- c(nnn, "x")
    fact <- factor(nnn, levels = nnn)
    ablabs <- character(n)
    ablabs[1] <- "alpha[1]"
    for (i in 1:(n-1)) ablabs[(i+1)] <- paste("beta[", i, "]", sep="")
    if (rel) for (i in 1:nrow(m)) {
        m[i,] <- m[i,] / m[i,ncol(m)]
        mp[i,] <- mp[i,] / mp[i,ncol(mp)]}
    hh <- c(m[1, 1], m[1,])
    hhh <- c(mp[1, 1], mp[1,])
    mh <- max(max(apply(m,1,max), apply(mp,1,max))*1.05, ymax)
    plot(fact, rep(-1000,length(fact)), type="n", 
        ylim=c(0,mh), xlim=c(-1, nrow(m)),
#        main=main, 
        xlab=xlab, ylab=ylab, ...)
    if (length(rownames(m)) == 1)
        fact <- factor(rownames(m), levels = rownames(m))
    for (i in 1:n) {
        for (j in 1:nrow(m)){
            if (!is.null(x$exp)) {
                if (plin[j,i]) {
                    lty <- 1 
                    lwd <- 5
                    col <- "black"
                    } else {
                    lty <- 1
                    lwd <- 5
                    col <- "grey"
                    }}
            if (i != n) arrows(j, m[j,i], j, m[j,(i+1)], code=0, lty=lty, lwd=lwd, col=col)
            if (j == 1) {
                arrows(-0.75, m[j,i], -1, m[j,i], code=0)
                if (!is.null(x$exp))
                    arrows(0.25, mp[j,i], 0, mp[j,i], code=0)
                legposo <- mean(hh[c(i, i+1)])
                legpose <- mean(hhh[c(i, i+1)])
                         if (i == 1) {
                            text(-0.5, m[j,i], expression(alpha[1]))
                            text(-0.5, mh, "Obs")
                            if (!is.null(x$exp)) text(0.5, mp[j,i], expression(alpha[1]))
                            } else {
                            text(-0.5, legposo, substitute(beta[z], list(z=i-1)))
                            if (!is.null(x$exp)) text(0.5, legpose, substitute(beta[z], list(z=i-1)))
                            }
                 }
            }
        if (attr(x, "design") == "twoway") {
            if (i != n) lines(fact, m[,i],type="l")
            if (i != n) lines(fact, mp[,i],type="l",lty=2)
            } else {
            lines(fact, m[,i],type="l")
            lines(fact, mp[,i],type="l",lty=2)}
        }
        for (i in n:1) {
        for (j in 1:nrow(m)){
            if (!is.null(x$exp)) {
                if (pdot[j,i]) {
                    pch <- 21
                    col2 <- "white"
                    bg <- "black"
                    } else {
                    pch <- 21
                    col2 <- "black"
                    bg <- "white"}
                }
            points(fact[j], m[j,i],type="p", pch=pch, bg=bg, col=col2, cex=1.2)
            }}

        if (!is.null(x$exp)){
            ptext <- paste("p <", attr(x, "crit"))
            text(0.5, mh, "Exp")
            lll <- strwidth(ptext)
            legend(p.legend, lty=1, lwd=lwd-2, col=c("black","grey"),
                legend=c(ptext, "NS"), xjust=1, text.width=lll, bty="n")
            legend(p.legend, pch=c(19, 19), col=c("black","white"),
                legend=c("", ""),bty="n", xjust=1, text.width=lll*1.1)
            legend(p.legend, pch=c(21, 21), col=c("white","black"),
                legend=c("", ""),bty="n", xjust=1, text.width=lll*1.1)}
}

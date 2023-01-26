`msoplot` <-
    function (x, alpha = 0.05, explained = FALSE, ylim = NULL,
    legend = "topleft", ...)
{
    if (is.data.frame(x$vario)) {
        vario <- x$vario
        hasSig <- is.numeric(x$vario$CA.signif)
        z <- qnorm(alpha/2)
        if (is.numeric(vario$CA.signif)) {
            vario <- vario[, -ncol(vario)]
        }
        ymax <- max(vario[, -1:-3], na.rm = TRUE)
        b <- ncol(vario) - 3
        label <- c("", "", "", "Total variance", "Explained plus residual",
                   "Residual variance", "Explained variance", "Conditioned variance")
        ci.lab <- "C.I. for total variance"
        sign.lab <- if(hasSig) "Sign. autocorrelation" else NULL
        if (is.numeric(x$CCA$rank)) {
            if (!explained)
                b <- b - 1
            if (is.numeric(x$vario$se))
                b <- b - 1
            figmat <- cbind(vario$All + z * vario$se,
                            vario$All - z * vario$se,
                            vario$Sum,
                            vario[, 6:(b + 3)])
            matplot(vario$Dist, cbind(0,figmat), type = "n",
                    xlab = "Distance", ylab = "Variance",
                    ylim = ylim, ...)
            lines(vario$Dist, vario$All + z * vario$se, lty = 1, ...)
            lines(vario$Dist, vario$All - z * vario$se, lty = 1, ...)
            lines(vario$Dist, vario$Sum, type = "b", lty = 2,
                  pch = 3, ...)
            ## Legend
            legend(legend,
                   legend=c(label[c(2,3:b)+3], ci.lab, sign.lab),
                   lty=c(c(1,2,1,1,1)[2:b], 1, if(hasSig) NA),
                   pch=c(3, (6:(b+3))-6, NA, if(hasSig) 15)
                   )
            matlines(vario$Dist, figmat[,-c(1:3)], type = "b", lty = 1,
                     pch = 6:(b+3)-6, ...)
            text(x = c(vario$Dist), y = par("usr")[3], pos = 3,
                 label = c(vario$n), cex = 0.8, ...)
            abline(v = max(x$H/2), lty = 3, ...)
        }
        else {
            if (is.null(ylim))
                ylim <- c(0, ymax)
            plot(vario$Dist, vario$All, type = "b", lty = 1,
                 pch = 0, xlab = "Distance", ylab = "Variance",
                 ylim = ylim, ...)
            abline(h = x$tot.chi, lty = 5, ...)
            text(x = c(vario$Dist), y = par("usr")[3], pos = 3,
                 label = c(vario$n), cex = 0.8)
            abline(v = max(x$H)/2, lty = 3, ...)
            legend(legend,
                   legend=c("Total variance","Global variance estimate",
                     if(hasSig) "Sign. autocorrelation"),
                   lty=c(1,5, if(hasSig) NA),
                   pch = if(hasSig) c(NA,NA,15) else NULL)
        }
    }
    if (hasSig) {
        a <- c(1:nrow(x$vario))[x$vario$CA.signif <
                                     alpha]
        points(vario$Dist[a], x$vario$CA[a], pch = 15, ...)
        if (is.numeric(x$CCA$rank)) {
            inflation <- 1 - weighted.mean(x$vario$CA, x$vario$n)/
                weighted.mean(x$vario$CA[-a],
                              x$vario$n[-a])
            cat("Error variance of regression model underestimated by",
                round(inflation * 100, 1), "percent", "\n")
        }
    }
    invisible()
}


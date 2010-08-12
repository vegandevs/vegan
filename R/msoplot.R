`msoplot` <-
    function (x, alpha = 0.05, explained = FALSE, ...) 
{
    object.cca <- x
    if (is.data.frame(object.cca$vario)) {
        object <- object.cca
        vario <- object$vario
        grain <- object$grain
        hasSig <- is.numeric(object$vario$CA.signif)
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
        ## You should not change par, or at least you must put
        ## back the old values when exiting:
        ## op <- par(omi = c(0.5, 0.5, 0, 0))
        ## on.exit(par(op))
        ##par(omi = c(0.5, 0.5, 0, 0))
        if (is.numeric(object$CCA$rank)) {
            if (!explained) 
                b <- b - 1
            if (is.numeric(object$vario$se)) 
                b <- b - 1
            plot(vario$Dist, vario$All, type = "n", lty = 1, 
                 pch = 3, xlab = "Distance", ylab = "Variance", 
                 ylim = c(0, ymax), cex.lab = 1.2, ...)
            lines(vario$Dist, vario$All + z * vario$se, lty = 1, ...)
            lines(vario$Dist, vario$All - z * vario$se, lty = 1, ...)
            lines(vario$Dist, vario$Sum, type = "b", lty = 2, 
                  pch = 3, ...)
            ## Legend
            legend("topleft", c(label[c(2,3:b)+3], ci.lab, sign.lab),
                   lty=c(c(1,2,1,1,1)[2:b], 1, if(hasSig) NA),
                   pch=c(3, (6:(b+3))-6, NA, if(hasSig) 15)
                   )
            for (i in 6:(b + 3)) {
                lines(vario$Dist, vario[, i], type = "b", lty = 1, 
                      pch = i - 6, ...)
            }
            text(x = c(vario$Dist), y = rep(0, length(vario$Dist)), 
                 label = c(vario$n), cex = 0.8, ...)
            lines(x = rep(max(object$H)/2, 2), y = c(-10, ymax + 
                                               10), lty = 3, ...)
        }
        else {
            plot(vario$Dist, vario$All, type = "b", lty = 1, 
                 pch = 0, xlab = "Distance", ylab = "Variance", 
                 ylim = c(0, ymax), cex.lab = 1.2, ...)
            lines(c(0, 10), rep(object$tot.chi, 2), lty = 5, ...)
            text(x = c(vario$Dist), y = rep(0, length(vario$Dist)), 
                 label = c(vario$n), cex = 0.8)
            lines(x = rep(max(object$H)/2, 2), y = c(-10, ymax + 
                                               10), lty = 3, ...)
            legend("topleft",
                   c("Total variance","Global variance estimate",
                     if(hasSig) "Sign. autocorrelation"),
                   lty=c(1,5, if(hasSig) NA),
                   pch = if(hasSig) c(NA,NA,15) else NULL)
        }
    }
    if (hasSig) {
        a <- c(1:nrow(object$vario))[object$vario$CA.signif < 
                                     alpha]
        points(vario$Dist[a], object$vario$CA[a], pch = 15, ...)
        if (is.numeric(object$CCA$rank)) {
            inflation <- 1 - weighted.mean(object$vario$CA, object$vario$n)/
                weighted.mean(object$vario$CA[-a], 
                              object$vario$n[-a])
            cat("Error variance of regression model underestimated by", 
                round(inflation * 100, 1), "percent", "\n")
        }
    }
    invisible()
}


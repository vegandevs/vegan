`plot.betadisper` <- function(x, axes = c(1,2), cex = 0.7, pch = seq_len(ng),
                              col = NULL, lty = "solid", lwd = 1, hull = TRUE,
                              ellipse = FALSE, conf,
                              segments = TRUE, seg.col = "grey",
                              seg.lty = lty, seg.lwd = lwd,
                              label = TRUE, label.cex = 1,
                              ylab, xlab, main, sub, ...)
{
    localAxis <- function(..., col, bg, pch, cex, lty, lwd) axis(...)
    localBox <- function(..., col, bg, pch, cex, lty, lwd) box(...)
    localTitle <- function(..., col, bg, pch, cex, lty, lwd) title(...)
    Ellipse <- function(scrs, centres, conf, col, lty, lwd, ...) {
        mat <- cov.wt(scrs, center = centres)
        if (mat$n.obs == 1)
            mat$cov[] <- 0
        xy <- if (mat$n.obs > 1) {
                  veganCovEllipse(mat$cov, mat$center, conf)
        } else {
            scrs
        }
        ordiArgAbsorber(xy, FUN = lines, col = col, lty = lty, lwd = lwd, ...)
    }
    if(missing(main))
        main <- deparse(substitute(x))
    if(missing(sub))
        sub <- paste("method = \"", attr(x, "method"), "\"", sep = "")
    if(missing(xlab))
        xlab <- paste("PCoA", axes[1])
    if(missing(ylab))
        ylab <- paste("PCoA", axes[2])
    t <- if (missing(conf)) {
        1
    } else {
        sqrt(qchisq(conf, df = 2))
    }
    g <- scores(x, choices = axes)
    ng <- length(levels(x$group))
    lev <- levels(x$group)
    ## sort out colour vector if none supplied
    if (is.null(col)) {
        col <- palette()
    }
    col <- rep_len(col, ng)        # make sure there are enough colors
    seg.col <- rep_len(seg.col, ng)     # ditto for segments
    plot(g$sites, asp = 1, type = "n", axes = FALSE, ann = FALSE, ...)
    ## if more than 1 group level
    if(is.matrix(g$centroids)) {
        for(i in seq_along(lev)) {
            curlev <- lev[i]
            take <- x$group == curlev
            j <- which(lev == curlev)
            if (segments) {
                segments(g$centroids[j, 1L], g$centroids[j, 2L],
                         g$sites[take, 1L],
                         g$sites[take, 2L], col = seg.col[i], lty = seg.lty,
                         lwd = seg.lwd)
            }
            if(hull) {
                ch <- chull(g$sites[take,, drop=FALSE])
                ch <- c(ch, ch[1])
                lines(x$vectors[take, axes, drop=FALSE][ch, ],
                      col = col[i], lty = lty, lwd = lwd, ...)
            }
            if (ellipse) {
                Ellipse(g$sites[take, , drop = FALSE],
                        centres = g$centroids[j, ],
                        conf = t,
                        col = col[i], lty = lty, lwd = lwd, ...)
            }
            points(g$centroids[j, , drop = FALSE], pch = 16, cex = 1,
                   col = col[i], ...)
        }
    } else {
        ## single group
        if (segments) {
            segments(g$centroids[, 1L], g$centroids[, 2L],
                     g$sites[, 1L], g$sites[, 2L], col = seg.col,
                     lty = seg.lty, ...)
        }
        if(hull) {
            ch <- chull(g$sites)
            ch <- c(ch, ch[1])
            lines(x$vectors[, axes][ch, ], col = col[1L], lty = lty,
                  lwd = lwd, ...)
        }
        if (ellipse) {
                Ellipse(g$sites,
                        centres = g$centroids,
                        conf = t,
                        col = col[1L], lty = lty, lwd = lwd,...)
        }
        points(g$centroids[, 1L], g$centroids[, 2L],
               pch = 16, cex = 1, col = col[1L], ...)
    }
    points(g$sites, pch = pch[x$group], cex = cex, col = col[x$group], ...)
    if (label) {
        ordilabel(x, display = "centroids", choices = axes, cex = label.cex,
            col = col)
    }
    localTitle(main = main, xlab = xlab, ylab = ylab, sub = sub, ...)
    localAxis(1, ...)
    localAxis(2, ...)
    localBox(...)
    class(g) <- "ordiplot"
    invisible(g)
}

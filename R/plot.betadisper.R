`plot.betadisper` <- function(x, axes = c(1,2), cex = 0.7, pch = seq_len(ng),
                              col = NULL, hull = TRUE,
                              ellipse = FALSE, ellipse.type = c("sd","se"),
                              ellipse.conf = NULL,
                              ylab, xlab, main, sub, ...)
{
    localAxis <- function(..., col, bg, pch, cex, lty, lwd) axis(...)
    localBox <- function(..., col, bg, pch, cex, lty, lwd) box(...)
    localTitle <- function(..., col, bg, pch, cex, lty, lwd) title(...)
    Ellipse <- function(scrs, type, conf, col, lty, lwd) {
        mat <- cov.wt(scrs)
        if (mat$n.obs == 1)
            mat$cov[] <- 0
        if (type == "se") {
            mat$cov <- mat$cov / mat$n.obs
        }
        xy <- if (mat$n.obs > 1) {
                  veganCovEllipse(mat$cov, mat$center, conf)
        } else {
            scrs
        }
        ordiArgAbsorber(xy, FUN = lines, col = col, lty = lty, lwd = lwd)
    }
    if(missing(main))
        main <- deparse(substitute(x))
    if(missing(sub))
        sub <- paste("method = \"", attr(x, "method"), "\"", sep = "")
    if(missing(xlab))
        xlab <- paste("PCoA", axes[1])
    if(missing(ylab))
        ylab <- paste("PCoA", axes[2])
    ellipse.type <- match.arg(ellipse.type)
    conf <- if (missing(ellipse.conf)) {
        1
    } else {
        sqrt(qchisq(conf, 2))
    }
    g <- scores(x, choices = axes)
    ng <- length(levels(x$group))
    lev <- levels(x$group)
    ## sort out colour vector if none supplied
    if (is.null(col)) {
        col <- palette()
    }
    col <- rep_len(col, ng)  # make sure there are enough colors
    plot(g$sites, asp = 1, type = "n", axes = FALSE, ann = FALSE, ...)
    ## if more than 1 group level
    if(is.matrix(g$centroids)) {
        for(i in seq_along(lev)) {
            curlev <- lev[i]
            take <- x$group == curlev
            j <- which(lev == curlev)
            segments(g$centroids[j, 1L], g$centroids[j, 2L],
                     g$sites[take, 1L],
                     g$sites[take, 2L], col = "blue", ...)
            if(hull) {
                ch <- chull(g$sites[take, ])
                ch <- c(ch, ch[1])
                lines(x$vectors[take, axes][ch, ],
                      col = col[i], lty = "dashed", ...)
            }
            if (ellipse) {
                Ellipse(g$sites[take, , drop = FALSE], type = ellipse.type,
                        conf = conf, col = col[i], lty = 1, lwd = 1)
            }
        }
        points(g$centroids, pch = 16, cex = 1, col = "red", ...)
    } else {
        ## single group
        segments(g$centroids[1L], g$centroids[2L],
                 g$sites[, 1L], g$sites[, 2L], col = "blue", ...)
        if(hull) {
            ch <- chull(g$sites)
            ch <- c(ch, ch[1])
            lines(x$vectors[, axes][ch, ],
                  col = col[1L], lty = "dashed", ...)
        }
        if (ellipse) {
                Ellipse(g$sites, type = ellipse.type,
                        conf = conf, col = col[1L], lty = 1, lwd = 1)
        }
        points(g$centroids[1L], g$centroids[1L],
               pch = 16, cex = 1, col = "red", ...)
    }
    points(g$sites, pch = pch[x$group], cex = cex, col = col[x$group], ...)
    localTitle(main = main, xlab = xlab, ylab = ylab, sub = sub, ...)
    localAxis(1, ...)
    localAxis(2, ...)
    localBox(...)
    class(g) <- "ordiplot"
    invisible(g)
}

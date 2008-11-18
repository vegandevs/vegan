`plot.cascadeKM` <-
    function (x, min.g, max.g, grpmts.plot = TRUE, sortg = FALSE, 
              gridcol = NA, ...) 
{
    wrapres <- x
    number <- (as.numeric(gsub(" groups", "", colnames(wrapres$results))))
    if (missing(min.g)) 
        min.g <- min(number)
    if (missing(max.g)) 
        max.g <- max(number)
    if (min.g < 2) 
        min.g = 2
    c.min <- which(number == min.g)
    c.max <- which(number == max.g)
    if (length(c.min) == 0) {
        stop("min.g value given has not been calculated by 'cascadeKM'\n")
    }
    if (length(c.max) == 0) {
        stop("max.g value given has not been calculated by 'cascadeKM'\n")
    }
    x <- wrapres$partition[, c.min:c.max]
    w <- wrapres$results[2, c.min:c.max]
    criterion <- wrapres$criterion
    x <- pregraphKM(x)
    if (sortg) {
        x <- orderingKM(x)
    }
    else {
    }
    main = (paste("K-means partitions comparison"))
    xlab = ("Number of groups in each partition")
    ylab = ("Objects")
    nc = ncol(x)
    colo <- (rainbow(max.g + 1))
    if (grpmts.plot) {
        def.par <- par(no.readonly = TRUE)
        nf <- layout(matrix(c(1, 2), 1, 2), widths = c(3, 1), 
                     TRUE)
        par(mar = c(5, 5, 5, 1), bg = "white", col = "black")
        image(1:nrow(x), 1:nc, x, col = colo, yaxt = "n", frame.plot = TRUE, 
              main = main, xlab = ylab, ylab = xlab, bg = NA)
        grid(nx = nrow(x), ny = max.g - min.g + 1, col = gridcol)
        box()
        axis(2, seq(min.g - min.g + 1, max.g - min.g + 1, by = 1), 
             lab = seq(min.g, max.g, by = 1))
        axis(1)
        par(mar = c(5, 2, 5, 1))
        par(bg = "white", fg = "black", col = "black")
        plot(y = min.g:max.g, x = w[1:nc], type = "b", main = paste(criterion, 
                                                       "\ncriterion", sep = ""), ylab = "", ylim = c(min.g - 
                                                                                            0.5, max.g + 0.5), yaxs = "i", yaxt = "n", xlab = "Values")
        grid(nx = NULL, ny = max.g - min.g + 1, col = gridcol)
        box()
        axis(2, seq(min.g, max.g, by = 1), lab = seq(min.g, max.g, 
                                           by = 1), col.axis = "black")
        axis(1)
        maxx = which.max(w[])
        minx = which.min(w[])
        tops <- which(w[c(2:nc)] - w[c(1:(nc - 1))] > 0) + 1
        maxx.o <- NA
        if (length(tops) != 0) {
            if (length(which(tops > maxx)) != 0) 
                maxx.o <- tops[which(tops > maxx)]
        }
        tops <- which(w[c(2:nc)] - w[c(1:(nc - 1))] < 0) + 1
        minx.o <- NA
        if (length(tops) != 0) {
            if (length(which(tops > minx)) != 0) 
                minx.o <- tops[which(tops > minx)]
        }
        if (tolower(criterion) == "calinski") {
            if (!is.na(maxx.o[1])) 
                points(y = maxx.o + min.g - 1, x = w[maxx.o], 
                       col = "orange", pch = 19)
            points(y = maxx + min.g - 1, x = w[maxx], col = "red", 
                   pch = 19)
        }
        else if (tolower(criterion) == "likelihood") {
            if (!is.na(maxx.o[1])) {
                points(y = maxx.o + min.g - 1, x = w[maxx.o], 
                       col = "orange", pch = 19)
            }
            points(y = maxx + min.g - 1, x = w[maxx], col = "red", 
                   pch = 19)
        }
        else if (tolower(criterion) == "ssi") {
            if (!is.na(maxx.o[1])) {
                points(y = maxx.o + min.g - 1, x = w[maxx.o], 
                       col = "orange", pch = 19)
            }
            points(y = maxx + min.g - 1, x = w[maxx], col = "red", 
                   pch = 19)
        }
        else {
            cat("When using the", criterion, "criterion, no red marker is", 
                "used to indicate the best value.\n")
        }
        par(def.par)
    }
    else {
        tops <- which(w[c(2:nc)] - w[c(1:(nc - 1))] > 0) + 1
        if (length(tops) != 0) {
            maxx <- which.max(w[c(2:nc)] - w[c(1:nc - 1)]) + 
                1
        }
        else {
            maxx <- which.max(w[])
            tops = 1
        }
    }
    out <- list(x = x, best.grps = maxx)
    if (grpmts.plot)
        invisible(out)
    else
        out
}

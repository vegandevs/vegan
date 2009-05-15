`taxa2dist` <-
    function (x, varstep = FALSE, check = TRUE, labels) 
{
    rich <- apply(x, 2, function(taxa) length(unique(taxa)))
    S <- nrow(x)
    if (check) {
        keep <- rich < S & rich > 1
        rich <- rich[keep]
        x <- x[, keep]
    }
    i <- rev(order(rich))
    x <- x[, i]
    rich <- rich[i]
    if (varstep) {
        add <- -diff(c(nrow(x), rich, 1))
        add <- add/c(S, rich)
        add <- add/sum(add) * 100
    }
    else {
        add <- rep(100/(ncol(x) + check), ncol(x) + check)
    }
    if (!is.null(names(add))) 
        names(add) <- c("Base", names(add)[-length(add)])
    if (!check) 
        add <- c(0, add)
    out <- matrix(add[1], nrow(x), nrow(x))
    for (i in 1:ncol(x)) {
        out <- out + add[i + 1] * outer(x[, i], x[, i], "!=")
    }
    out <- as.dist(out)
    attr(out, "method") <- "taxa2dist"
    attr(out, "steps") <- add
    if (missing(labels)) {
        attr(out, "Labels") <- rownames(x)
    } else {
        if (length(labels) != nrow(x))
            warning("Labels are wrong: needed ", nrow(x), " got ", length(labels))
        attr(out, "Labels") <- as.character(labels)
    }
    if (!check && any(out <= 0))
        warning("you used 'check=FALSE' and some distances are zero -- was this intended?")
    out
}

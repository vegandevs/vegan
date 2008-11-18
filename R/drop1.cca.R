`drop1.cca` <-
    function(object, scope, test = c("none", "permutation"),
             pstep = 100, perm.max = 200, ...)
{
    test <- match.arg(test)
    out <- NextMethod("drop1", object, test="none", ...)
    cl <- class(out)
    if (test == "permutation") {
        rn <- rownames(out)[-1]
        if (missing(scope))
            scope <- rn
        adds <- anova(object, by = "margin", step = pstep,
                      perm.max = perm.max, scope = scope, ...)
        nr <- nrow(adds)
        out <- cbind(out, rbind(NA, adds[rn,3:5]))
        class(out) <- cl
    }
    out
}

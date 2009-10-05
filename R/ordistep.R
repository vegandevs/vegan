`ordistep` <-
    function(object, scope, Pin = 0.05, Pout = 0.1,
             pstep = 100, perm.max = 1000, steps=50, trace = TRUE, ...)
{
    if (!inherits(object, "cca"))
        stop("function can be only used to 'cca' and related objects")
    if (missing(scope))
        stop("absolutely needs scope")
    mod <- eval.parent(update(object, . ~  1))
    for (i in 1:steps){
        change <- NULL
        ## Consider dropping
        if (length(drop.scope(mod))) {
            aod <- drop1(mod, test="perm", pstep = pstep, perm.max = perm.max,
                         alpha = Pout)
            aod <- aod[-1,]
            o <- order(-aod[,5], aod[,4], aod[,2])
            aod <- aod[o,]
            rownames(aod) <- paste("-", rownames(aod), sep = " ")
            if (trace) {
                cat("\n")
                print(aod)
            }
            if (aod[1,5] > Pout) {
                change <- rownames(aod)[1]
                mod <- eval.parent(update(mod, paste("~  .", change)))
                if (trace) {
                    cat("\n")
                    print(mod$call)
                }
            }
        }
        ## Consider adding
        if (length(add.scope(formula(mod), scope))) {
            aod <- add1(mod, scope = scope, test = "perm", pstep = pstep,
                        perm.max = perm.max, alpha = Pin)
            aod <- aod[-1,]
            o <- order(aod[,5], aod[,4], aod[,2])
            aod <- aod[o,]
            rownames(aod) <- paste("+", rownames(aod), sep = " ")
            if (trace) {
                cat("\n")
                print(aod)
            }
            if (aod[1,5] <= Pin) {
                change <- rownames(aod)[1]
                mod <- eval.parent(update(mod, paste( "~  .",change)))
                if (trace) {
                    cat("\n")
                    print(mod$call)
                }
            }
        }
        ## No drop, no add: done
        if (is.null(change))
            break
    }
    cat("\n")
    mod
}

`ordistep` <-
    function(object, scope, direction =c("both", "backward", "forward"),
             Pin = 0.05, Pout = 0.1, permutations = how(nperm = 199),
             steps=50, trace = TRUE, ...)
{
    if (!inherits(object, "cca"))
        stop("function can be only used with 'cca' and related objects")
    if (is.null(object$terms))
        stop("ordination model must be fitted using formula")
    ## handling 'direction' and 'scope' directly copied from
    ## stats::step()
    md <- missing(direction)
    direction <- match.arg(direction)
    backward <- direction == "both" | direction == "backward"
    forward <- direction == "both" | direction == "forward"
    ffac <- attr(terms(object), "factors")
    if (missing(scope)) {
        fdrop <- numeric(0)
        fadd <- ffac
        if (md)
            forward <- FALSE
    }
    else {
        if (is.list(scope) && (!is.null(scope$lower) || !is.null(scope$upper))) {
            fdrop <- if (!is.null(fdrop <- scope$lower))
                attr(terms(update.formula(object, fdrop)), "factors")
            else numeric(0)
            fadd <- if (!is.null(fadd <- scope$upper))
                attr(terms(update.formula(object, fadd)), "factors")
        }
        else {
            fadd <- if (!is.null(fadd <- scope))
                        attr(terms(update.formula(object, scope)), "factors")
            if (forward)
                fdrop <- attr(terms(object), "factor")
            else
                fdrop <- numeric(0L)
        }
    }
    scope <- factor.scope(ffac, list(add = fadd, drop = fdrop))
    ## 'anotab' collects the changes into 'anova' object in the output
    anotab <- NULL
    if (trace) {
        cat("\n")
        cat(pasteCall(formula(object), prefix = "Start:"))
    }
    for (i in 1:steps){
        change <- NULL
        ## Consider dropping
        if (backward && length(scope$drop)) {
            aod <- drop1(object, scope = scope$drop, test="perm",
                         permutations = permutations,
                         alpha = Pout, trace = trace, ...)
            aod <- aod[-1,]
            o <- order(-aod[,4], aod[,2])
            aod <- aod[o,]
            rownames(aod) <- paste("-", rownames(aod), sep = " ")
            if (trace) {
                cat("\n")
                print(aod)
            }
            if (is.na(aod[1,4]) || aod[1,4] > Pout) {
                anotab <- rbind(anotab, aod[1,])
                change <- rownames(aod)[1]
                object <- eval.parent(update(object, paste("~  .", change)))
                scope <- factor.scope(attr(terms(object), "factors"),
                                      list(add = fadd, drop = fdrop))
                if (trace) {
                    cat("\n")
                    cat(pasteCall(formula(object), prefix = "Step:"))
                }
            }
        }
        ## Consider adding
        if (forward && length(scope$add)) {
            aod <- add1(object, scope = scope$add, test = "perm",
                        permutations = permutations,
                        alpha = Pin, trace = trace, ...)
            aod <- aod[-1,]
            o <- order(aod[,4], aod[,2])
            aod <- aod[o,]
            rownames(aod) <- paste("+", rownames(aod), sep = " ")
            if (trace) {
                cat("\n")
                print(aod)
            }
            if (!is.na(aod[1,4]) && aod[1,4] <= Pin) {
                anotab <- rbind(anotab, aod[1,])
                change <- rownames(aod)[1]
                object <- eval.parent(update(object, paste( "~  .",change)))
                scope <- factor.scope(attr(terms(object), "factors"),
                                      list(add = fadd, drop = fdrop))
                if (trace) {
                    cat("\n")
                    cat(pasteCall(formula(object), prefix="Step:"))
                }
            }
        }
        ## No drop, no add: done
        if (is.null(change))
            break
    }
    cat("\n")
    object$anova <- anotab
    object
}

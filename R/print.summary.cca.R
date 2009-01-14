`print.summary.cca` <-
    function (x, digits = x$digits, head=NA, tail=head, ...) 
{
    hcat <- function(x, head=head, tail=tail, ...) {
        if (!is.na(head) && !is.na(tail) && head+tail+4 < nrow(x))
            x <- rbind(head(x, n=head), "...." = NA, tail(x, n=tail))
        print(x, na.print = "", ...)
    }
    cat("\nCall:\n")
    cat(deparse(x$call), "\n")
    statnam <- if (x$method == "cca") 
        "averages"
    else "sums"
    cat("\nPartitioning of ", x$inertia, ":\n", sep = "")
    out <- c(Total = x$tot.chi, Conditioned = x$partial.chi, 
             Constrained = x$constr.chi, Unconstrained = x$unconst.chi)
    out <- cbind(Inertia = out, Proportion = out/out[1])
    print(out, digits = digits, ...)
    cat("\nEigenvalues, and their contribution to the", x$inertia, 
        "\n")
    if (!is.null(x$partial.chi)) {
        cat("after removing the contribution of conditiniong variables\n")
    }
    cat("\n")
    out <- rbind("Eig.value" = c(x$ev.con, x$ev.uncon),
                 "Accounted" = c(x$ev.con.account, x$ev.uncon.account))
    print(out, digits = digits, ...)
    if (!is.null(x$cca.acc)) {
        cat("\nAccumulated constrained eigenvalues\n")
        print(x$cca.acc, digits = digits, ...)
    }
    cat("\nScaling", x$scaling, "for species and site scores\n")
    if (abs(x$scaling) == 2) {
        ev.ent <- "Species"
        other.ent <- "Sites"
    }
    else if (abs(x$scaling) == 1) {
        ev.ent <- "Sites"
        other.ent <- "Species"
    }
    else if (abs(x$scaling) == 3) {
        ev.ent <- "Both sites and species"
        other.ent <- NULL
    }
    if (x$scaling) {
        cat("*", ev.ent, "are scaled proportional to eigenvalues\n")
        if (!is.null(other.ent)) 
            cat("*", other.ent, "are unscaled: weighted dispersion equal")
        cat(" on all dimensions\n")
    }
    if (!x$scaling) {
        cat("* Both are 'unscaled' or as they are in the result\n")
    }
    if (x$scaling < 0) {
        if (x$method == "cca") 
            cat("* Hill scaling performed on both scores\n")
        else 
            cat("* Species scores divided by species standard deviations\n")
        cat("  so that they no longer are biplot scores\n")
    }
    if (x$method != "cca") {
        cat("* General scaling constant of scores: ",
            attr(x, "const"), "\n")
    }
    if (!is.null(x$species)) {
        cat("\n\nSpecies scores\n\n")
        hcat(x$species, head=head, tail=tail, digits = digits, ...)
    }
    if (!is.null(x$sites)) {
        cat("\n\nSite scores (weighted", statnam, "of species scores)\n\n")
        hcat(x$sites, head=head, tail=tail, digits = digits, ...)
    }
    if (!is.null(x$constraints)) {
        cat("\n\nSite constraints (linear combinations of constraining variables)\n\n")
        hcat(x$constraints, head=head, tail=tail, digits = digits, ...)
    }
    if (!is.null(x$biplot)) {
        cat("\n\nBiplot scores for constraining variables\n\n")
        print(x$biplot, digits = digits, ...)
    }
    if (!is.null(x$centroids) && !is.na(x$centroids[1])) {
        cat("\n\nCentroids for factor constraints\n\n")
        print(x$centroids, digits = digits, ...)
    }
    cat("\n")
    invisible(x)
}

### Make a compact summary of permutations. This copies Gav Simpson's
### permute:::print.how, but only displays non-default choices in how().
`howHead` <- function(x, ...)
{
    ## print nothing is this not 'how'
    if (is.null(x) || !inherits(x, "how"))
        return()
    ## collect header
    head <- NULL
    ## blocks
    if (!is.null(getBlocks(x)))
        head <- paste0(head, paste("Blocks: ", x$blocks.name, "\n"))
    ## plots
    plotStr <- getStrata(x, which = "plots")
    if (!is.null(plotStr)) {
        plots <- getPlots(x)
        ptype <- getType(x, which = "plots")
        head <- paste0(head, paste0("Plots: ", plots$plots.name, ", "))
        head <- paste0(head, paste("plot permutation:", ptype))
        if(getMirror(x, which = "plots"))
            head <- paste(head, "mirrored")
        if (ptype == "grid") {
            nr <- getRow(x, which = "plots")
            nc <- getCol(x, which = "plots")
            head <- paste0(head, sprintf(ngettext(nr, " %d row", " %d rows"),
                                        nr))
            head <- paste0(head, sprintf(ngettext(nc, " %d column",
                                                 " %d columns"), nc))
        }
        head <- paste0(head, "\n")
    }
    ## the fine level (within plots if any)
    type <- getType(x, which = "within")
    head <- paste0(head, "Permutation: ", type)
    if (type %in% c("series", "grid")) {
        if(getMirror(x, which = "within"))
            head <- paste(head, "mirrored")
        if(getConstant(x))
            head <- paste0(head, " constant permutation within each Plot")
    }
    if (type == "grid") {
        nr <- getRow(x, which = "within")
        nc <- getCol(x, which = "within")
        head <- paste0(head, sprintf(ngettext(nr, " %d row", " %d rows"),
                                    nr))
        head <- paste0(head, sprintf(ngettext(nc, " %d column",
                                             " %d columns"), nc))
    }
    paste0(head, "\nNumber of permutations: ", getNperm(x),  "\n")
}

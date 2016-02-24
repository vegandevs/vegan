`vegandocs` <-
    function (doc = c("NEWS", "ONEWS", "FAQ-vegan",
              "intro-vegan", "diversity-vegan",
              "decision-vegan", "partitioning", "permutations")) 
{
    doc <- match.arg(doc)
    if (doc == "NEWS")
        news(package = "vegan")
    else if (doc %in% vignette(package="vegan")$results[, "Item"])
        vignette(doc, package = "vegan")
    else if (doc == "permutations")
        vignette(doc, package = "permute")
    else if (doc == "partitioning") { # pdf but not a vignette
        doc <- paste0(doc, ".pdf")
        doc <- file.path(system.file(package="vegan"), "doc", doc)
        if (.Platform$OS.type == "windows")
            shell.exec(doc)
        else system(paste(getOption("pdfviewer"), doc, "&"))
    }
    else # last resort
        file.show(system.file(package="vegan", doc))
}

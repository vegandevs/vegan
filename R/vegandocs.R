`vegandocs` <-
    function (doc = c("NEWS", "ONEWS", "FAQ-vegan",
              "intro-vegan", "diversity-vegan",
              "decision-vegan", "partitioning", "permutations"))
{
    doc <- match.arg(doc)
    if (doc == "NEWS") {
        .Defunct('news(package="vegan")')
    } else if (doc %in% vignette(package="vegan")$results[, "Item"]) {
        .Defunct('browseVignettes("vegan")')
    } else if (doc == "permutations") {
        .Defunct('browseVignettes("permute")')
    } else { # last resort
        .Defunct(gettextf('file.show(system.file(package="vegan", "%s"))', doc))
    }
}

`vegandocs` <-
    function (doc = c("NEWS", "ONEWS", "FAQ-vegan",
              "intro-vegan", "diversity-vegan",
              "decision-vegan", "partitioning", "permutations"))
{
    doc <- match.arg(doc)
    if (doc == "NEWS") {
        .Defunct('news(package="vegan")')
        news(package = "vegan")
    } else if (doc %in% vignette(package="vegan")$results[, "Item"]) {
        .Defunct('browseVignettes("vegan")')
        vignette(doc, package = "vegan")
    } else if (doc == "permutations") {
        .Defunct('browseVignettes("permute")')
        vignette(doc, package = "permute")
    } else  # last resort
        file.show(system.file(package="vegan", doc))
}

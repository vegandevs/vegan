`vegandocs` <-
    function (doc = c("NEWS", "ONEWS", "ChangeLog", "FAQ-vegan.pdf",
              "intro-vegan.pdf", "diversity-vegan.pdf",
              "decision-vegan.pdf", "partitioning.pdf", "permutations.pdf")) 
{
    doc <- match.arg(doc)
    if (length(grep(".pdf", doc)) > 0) {
        if (doc == "permutations.pdf")
            doc <- file.path(system.file(package="permute"), "doc", doc)
        else
            doc <- file.path(system.file(package="vegan"), "doc", doc)
        if (.Platform$OS.type == "windows")
            shell.exec(doc)
        else system(paste(getOption("pdfviewer"), doc, "&"))
    } else if (doc == "NEWS") {
        ## Try html
        helptype <- getOption("help_type")
        if (length(helptype) && helptype == "html") {
            browseURL(paste("file://",
                            system.file(package="vegan", "doc", "NEWS.html"),
                            sep=""))
        } else {
            file.show(Rd2txt(file.path(system.file(package="vegan"),
                                               "NEWS.Rd"), tempfile()))
        }
    } else {
        file.show(system.file(package="vegan", doc))
    } 
}

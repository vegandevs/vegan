###
### Editable Tcl/Tk plot for ordination
###
`orditkplot` <-
    function(x, display = "species", choices = 1:2, width, xlim, ylim,
             tcex=0.8, tcol, pch = 1, pcol, pbg, pcex = 0.7,
             labels,  ...)
{
    if (!capabilities("tcltk"))
        stop("Your R has no capability for Tcl/Tk")
    requireNamespace("tcltk") || stop("requires package tcltk")

############################
### Check and sanitize input
###########################

    ## Graphical parameters and constants, and save some for later plotting
    p <- par()
    sparnam <- c("bg","cex", "cex.axis","cex.lab","col", "col.axis", "col.lab",
                 "family", "fg", "font", "font.axis", "font.lab", "lheight",
                 "lwd", "mar", "mex", "mgp", "ps", "tcl", "las")
    ## Get par given in the command line and put them to p
    if (inherits(x, "orditkplot")) {
        dots <- x$par
        for (arg in names(x$args))
            assign(arg, unlist(x$args[arg]))
    } else {
        dots <- match.call(expand.dots = FALSE)$...
    }
    if (length(dots) > 0) {
        dots <- dots[names(dots) %in% sparnam]
        ## eval() or mar=c(4,4,1,1) will be a call, not numeric
        dots <- lapply(dots, function(x) if (is.call(x)) eval(x) else x)
        p <- check.options(new = dots, name.opt = "p",
                           envir = environment())
    }
    savepar <- p[sparnam]
    PPI <- 72                                         # Points per Inch
    p2p <- as.numeric(tcltk::tclvalue(tcltk::tcl("tk", "scaling"))) # Pixel per point
    DIAM <- 2.7                               # diam of plotting symbol
    ## Plotting symbol diam
    diam <- round(pcex * DIAM * p2p, 1)
    ## Sanitize colours
    sanecol <- function(x) {
        if (is.numeric(x))
            x <- palette()[x]
        x <- gsub("transparent", "", x)
        x[is.na(x)] <- ""
        x
    }
    p$bg <- sanecol(p$bg)
    p$fg <- sanecol(p$fg)
    p$col <- sanecol(p$col)
    p$col.axis <- sanecol(p$col.axis)
    p$col.lab <- sanecol(p$col.lab)
    ## Point and label colours
    if (missing(pcol))
        pcol <- p$col
    if (missing(pbg))
        pbg <- "transparent"
    if (missing(tcol))
        tcol <- p$col
    pcol <- sanecol(pcol)
    pbg <- sanecol(pbg)
    tcol <- sanecol(tcol)
    ## Define fonts
    idx <- match(p$family, c("","serif","sans","mono"))
    if (!is.na(idx))
        p$family <- c("Helvetica", "Times", "Helvetica", "Courier")[idx]
    saneslant <- function(x) {
        list("roman", "bold", "italic", c("bold", "italic"))[[x]]
    }
    ## fnt must be done later, since family, font and size can be
    ## vectors and slant can be of length 1 or 2
    ## fnt <- c(p$family, round(p$ps*p$cex*tcex), saneslant(p$font))
    labfam <- p$family
    labsize <- round(p$ps * p$cex * tcex)
    fnt.axis <- c(p$family, round(p$ps*p$cex.axis), saneslant(p$font.axis))
    fnt.lab <- c(p$family, round(p$ps*p$cex.lab), saneslant(p$font.lab))
    ## Imitate R plotting symbols pch
    SQ <- sqrt(2)     # Scaling factor for plot
    Point <- function(x, y, pch, col, fill, diam) {
        x <- round(x)
        y <- round(y)
        switch(as.character(pch),
               "0" = Point(x, y, 22, col, fill = "", diam),
               "1" = Point(x, y, 21, col, fill = "", diam),
               "2" = Point(x, y, 24, col, fill = "", diam),
               "3" = {tcltk::tkcreate(can, "line",
                               x, y+SQ*diam, x, y-SQ*diam, fill=col)
                      tcltk::tkcreate(can, "line",
                               x+SQ*diam, y, x-SQ*diam, y, fill=col)},
               "4" = {tcltk::tkcreate(can, "line",
                               x-diam, y-diam, x+diam, y+diam, fill=col)
                      tcltk::tkcreate(can, "line",
                               x-diam, y+diam, x+diam, y-diam, fill=col)},
               "5" = Point(x, y, 23, col, fill = "", diam),
               "6" = Point(x, y, 25, col, fill = "", diam),
               "7" = {Point(x, y, 4, col, fill, diam)
                      Point(x, y, 0, col, fill, diam)},
               "8" = {Point(x, y, 3, col, fill, diam)
                      Point(x, y, 4, col, fill, diam)},
               "9" = {Point(x, y, 3, col, fill, diam)
                      Point(x, y, 5, col, fill, diam)},
               "10" = {Point(x, y, 3, col, fill, diam/SQ)
                       Point(x, y, 1, col, fill, diam)},
               "11" = {Point(x, y, 2, col, fill, diam)
                       Point(x, y, 6, col, fill, diam)},
               "12" = {Point(x, y, 3, col, fill, diam/SQ)
                       Point(x, y, 0, col, fill, diam)},
               "13" = {Point(x, y, 4, col, fill, diam)
                       Point(x, y, 1, col, fill, diam)},
               "14" = {tcltk::tkcreate(can, "line", x-diam, y-diam, x, y+diam,
                                fill = col)
                       tcltk::tkcreate(can, "line", x+diam, y-diam, x, y+diam,
                                fill = col)
                       Point(x, y, 0, col, fill, diam)},
               "15" = Point(x, y, 22, col = col, fill = col, diam),
               "16" = Point(x, y, 21, col = col, fill = col, diam),
               "17" = Point(x, y, 24, col = col, fill = col, diam),
               "18" = Point(x, y, 23, col = col, fill = col, diam/SQ),
               "19" = Point(x, y, 21, col = col, fill = col, diam),
               "20" = Point(x, y, 21, col = col, fill = col, diam/2),
               "21" = tcltk::tkcreate(can, "oval", x-diam, y-diam,
               x+diam, y+diam, outline = col, fill = fill),
               "22" = tcltk::tkcreate(can, "rectangle", x-diam, y-diam,
               x+diam, y+diam, outline = col, fill = fill),
               "23" = tcltk::tkcreate(can, "polygon", x, y+SQ*diam,
               x+SQ*diam, y, x, y-SQ*diam, x-SQ*diam, y,
               outline = col, fill = fill),
               "24" = tcltk::tkcreate(can, "polygon", x, y-SQ*diam,
               x+sqrt(6)/2*diam, y+SQ/2*diam, x-sqrt(6)/2*diam, y+SQ/2*diam,
               outline = col, fill = fill),
               "25" = tcltk::tkcreate(can, "polygon", x, y+SQ*diam,
               x+sqrt(6)/2*diam, y-SQ/2*diam, x-sqrt(6)/2*diam, y-SQ/2*diam,
               outline = col, fill = fill),
               "o" = Point(x, y, 1, col, fill, diam),
               ## default: text with dummy location of the label
               {tcltk::tkcreate(can, "text",
                        x, y, text = as.character(pch), fill = col)
                Point(x, y, 21, col="", fill="", diam)}
                )
    }

############################
### Initialize Tcl/Tk Window
############################

    ## toplevel
    w <- tcltk::tktoplevel()
    tcltk::tktitle(w) <- deparse(match.call())
    ## Max dim of windows (depends on screen)
    YSCR <- as.numeric(tcltk::tkwinfo("screenheight", w)) - 150
    XSCR <- as.numeric(tcltk::tkwinfo("screenwidth", w))

################################
### Buttons and button functions
################################

    ## Buttons
    buts <- tcltk::tkframe(w)
    ## Copy current canvas to EPS using the standard Tcl/Tk utility
    cp2eps <- tcltk::tkbutton(buts, text="Copy to EPS",
                       command=function() tcltk::tkpostscript(can, x=0, y=0,
                       height=height, width=width,
                       file=tcltk::tkgetSaveFile(
                         filetypes="{{EPS file} {.eps}}",
                         defaultextension=".eps")))
    dismiss <- tcltk::tkbutton(buts, text="Dismiss",
                               command=function() tcltk::tkdestroy(w))
    ## Dump current plot to an "orditkplot" object (internally)
    ordDump <- function() {
        xy <- matrix(0, nrow=nrow(sco), ncol=2)
        rownames(xy) <- rownames(sco)
        colnames(xy) <- colnames(sco)
        for(nm in names(pola)) {
            xy[as.numeric(tcltk::tclvalue(id[[nm]])),] <- xy2usr(nm)
        }
        curdim <- round(c(width, height) /PPI/p2p, 2)
        ## Sanitize colours for R plot
        pbg[pbg == ""]  <- "transparent"
        pcol[pcol == ""] <- "transparent"
        ## Reduce vector args if all entries are constant
        argcollapse <- function(x)
            if (length(unique(x)) == 1) x[1] else x
        pch <- argcollapse(pch)
        pcol <- argcollapse(pcol)
        pbg <- argcollapse(pbg)
        tcol <- argcollapse(tcol)
        ## Save
        args <- list(tcex = tcex, tcol = tcol, pch = pch, pcol = pcol,
                     pbg = pbg, pcex = pcex, xlim = xlim, ylim = ylim)
        xy <- list(labels = xy, points = sco, par = savepar, args = args,
                   dim = curdim)
        class(xy) <- "orditkplot"
        xy
    }
    ## Button to dump "orditkplot" object to the R session
    pDump <- function() {
        xy <- ordDump()
        dumpVar <- tcltk::tclVar("")
        tt <- tcltk::tktoplevel()
        tcltk::tktitle(tt) <- "R Dump"
        entryDump <- tcltk::tkentry(tt, width=20, textvariable=dumpVar)
        tcltk::tkgrid(tcltk::tklabel(tt, text="Enter name for an R object"))
        tcltk::tkgrid(entryDump, pady="5m")
        isDone <- function() {
            dumpName <- tcltk::tclvalue(dumpVar)
            if (exists(dumpName, envir=.GlobalEnv)) {
                ok <- tcltk::tkmessageBox(message=paste(sQuote(dumpName),
                                   "exists.\nOK to overwrite?"),
                                   icon="warning", type="okcancel",
                                   default="ok")
                if(tcltk::tclvalue(ok) == "ok") {
                    assign(dumpName, xy, envir=.GlobalEnv)
                    tcltk::tkdestroy(tt)
                }
            }
            else {
                assign(dumpName, xy, envir=.GlobalEnv)
                tcltk::tkdestroy(tt)
            }
        }
        tcltk::tkbind(entryDump, "<Return>", isDone)
        tcltk::tkfocus(tt)
    }
    dump <- tcltk::tkbutton(buts, text="Dump to R", command=pDump)
    ## Button to write current "orditkplot" object to a graphical device
    devDump <- function() {
        xy <- ordDump()
        ftypes <- c("eps" = "{EPS File} {.eps}",
                    "pdf" = "{PDF File} {.pdf}",
                    "svg" = "{SVG File} {.svg}",
                    "png" = "{PNG File} {.png}",
                    "jpg" = "{JPEG File} {.jpg .jpeg}",
                    "bmp" = "{BMP File} {.bmp}",
                    "tiff"= "{TIFF File} {.tif .tiff}",
                    "fig" = "{XFig File} {.fig}")
        falt <- rep(TRUE, length(ftypes))
        names(falt) <- names(ftypes)
        if (!capabilities("png"))
            falt["png"] <- FALSE
        if (!capabilities("jpeg"))
            falt["jpg"] <- FALSE
        if (!capabilities("cairo"))
            falt["svg"] <- FALSE
        ## Should work also in R < 2.8.0 with no capabilities("tiff")
        if (!isTRUE(unname(capabilities("tiff"))))
            falt["tiff"] <- FALSE
        ftypes <- ftypes[falt]
        ## External Tcl/Tk in Windows seems to buggy with type
        ## extensions of the file name: the extension is not
        ## automatically appended, but defaultextension is interpreted
        ## wrongly so that its value is not used as extension but
        ## correct appending is done if defaultextension has any
        ## value. The following kluge is against Tcl/Tk documentation,
        ## and should be corrected if Tcl/Tk is fixed.
        if (.Platform$OS.type == "windows")
            fname <- tcltk::tkgetSaveFile(filetypes=ftypes,
                                          defaultextension = TRUE)
        else
            fname <- tcltk::tkgetSaveFile(filetypes=ftypes)
        if(tcltk::tclvalue(fname) == "")
            return(NULL)
        fname <- tcltk::tclvalue(fname)
        ftype <- unlist(strsplit(fname, "\\."))
        ftype <- ftype[length(ftype)]
        if (ftype == "jpeg")
            ftype <- "jpg"
        if (ftype == "tif")
            ftype <- "tiff"
        mess <- "is not a supported type: file not produced. Supported types are"
        if (!(ftype %in% names(ftypes))) {
            tcltk::tkmessageBox(message=paste(sQuote(ftype), mess, paste(names(ftypes),
                         collapse=", ")), icon="warning")
            return(NULL)
        }
        pixdim <- round(xy$dim*PPI*p2p)
        switch(ftype,
               eps = postscript(file=fname, width=xy$dim[1], height=xy$dim[2],
               paper="special", horizontal = FALSE),
               pdf = pdf(file=fname, width=xy$dim[1], height=xy$dim[2]),
               svg = svg(filename=fname, width=xy$dim[1], height=xy$dim[2]),
               png = png(filename=fname, width=pixdim[1], height=pixdim[2]),
               jpg = jpeg(filename=fname, width=pixdim[1], height=pixdim[2],
               quality = 100),
               tiff = tiff(filename=fname, width=pixdim[1], height=pixdim[2]),
               bmp = bmp(filename=fname, width=pixdim[1], height=pixdim[2]),
               fig = xfig(file=fname, width=xy$dim[1], height=xy$dim[2]))
        plot.orditkplot(xy)
        dev.off()
    }
    export <- tcltk::tkbutton(buts, text="Export plot", command=devDump)

##########
### Canvas
##########

    ## Make canvas
    sco <- try(scores(x, display=display, choices = choices, ...),
               silent = TRUE)
    if (inherits(sco, "try-error")) {
        tcltk::tkmessageBox(message=paste("No ordination scores were found in",
                     sQuote(deparse(substitute(x)))), icon="error")
        tcltk::tkdestroy(w)
        stop("argument x did not contain ordination scores")
    }
    if (!missing(labels))
        rownames(sco) <- labels
    ## Recycle graphical parameters in plots
    nr <- nrow(sco)
    pcol <- rep(pcol, length=nr)
    pbg <- rep(pbg, length=nr)
    pch <- rep(pch, length=nr)
    tcol <- rep(tcol, length=nr)
    diam <- rep(diam, length=nr)
    labfam <- rep(labfam, length=nr)
    labsize <- rep(labsize, length=nr)
    if (inherits(x, "ordipointlabel"))
        labfnt <- attr(x$labels, "font")
    else
        labfnt <- rep(p$font, length=nr)
    ## Select only items within xlim, ylim
    take <- rep(TRUE, nr)
    if (!missing(xlim))
        take <- take & sco[,1] >= xlim[1] & sco[,1] <= xlim[2]
    if (!missing(ylim))
        take <- take & sco[,2] >= ylim[1] & sco[,2] <= ylim[2]
    sco <- sco[take,, drop=FALSE]
    labs <- rownames(sco)
    pcol <- pcol[take]
    pbg <- pbg[take]
    tcol <- tcol[take]
    pch <- pch[take]
    diam <- diam[take]
    labfam <- labfam[take]
    labsize <- labsize[take]
    labfnt <- labfnt[take]
    ## Ranges and pretty values for axes
    if (missing(xlim))
        xlim <- range(sco[,1], na.rm = TRUE)
    if (missing(ylim))
        ylim <- range(sco[,2], na.rm = TRUE)
    xpretty <- pretty(xlim)
    ypretty <- pretty(ylim)
    ## Extend ranges by 4%
    xrange <- c(-0.04, 0.04) * diff(xlim) + xlim
    xpretty <- xpretty[xpretty >= xrange[1] & xpretty <= xrange[2]]
    yrange <- c(-0.04, 0.04) * diff(ylim) + ylim
    ypretty <- ypretty[ypretty >= yrange[1] & ypretty <= yrange[2]]
    ## Canvas like they were in the default devices when I last checked
    if (missing(width))
        width <- p$din[1]
    width <- width * PPI * p2p
    ## Margin row width also varies with platform and devices
    ## rpix <- (p$mai/p$mar * PPI * p2p)[1]
    rpix <- p$cra[2]
    mar <- round(p$mar * rpix)
    xusr <- width - mar[2] - mar[4]
    xincr <- xusr/diff(xrange)
    yincr <- xincr
    xy0 <- c(xrange[1], yrange[2]) # upper left corner
    ## Functions to translate scores to canvas coordinates and back
    usr2xy <- function(row) {
        x <- (row[1] - xy0[1]) * xincr + mar[2]
        y <- (xy0[2] - row[2]) * yincr + mar[3]
        c(x,y)
    }
    ## User coordinates of an item
    xy2usr <- function(item) {
        xy <- as.numeric(tcltk::tkcoords(can, item))
        x <- xy[1]
        y <- xy[2]
        x <- xrange[1] + (x - mar[2])/xincr
        y <- yrange[2] - (y - mar[3])/yincr
        c(x,y)
    }
    ## Canvas x or y to user coordinates
    x2usr <- function(xcan) {
        xrange[1] + (xcan - mar[2])/xincr
    }
    y2usr <- function(ycan) {
        yrange[2] - (ycan - mar[3])/yincr
    }
    ## Equal aspect ratio
    height <- round((diff(yrange)/diff(xrange)) * xusr)
    height <- height + mar[1] + mar[3]
    ## Canvas, finally
    can <- tcltk::tkcanvas(w, relief="sunken", width=width, height=min(height,YSCR),
                    scrollregion=c(0,0,width,height))
    if (p$bg != "")
        tcltk::tkconfigure(can, bg=p$bg)
    yscr <- tcltk::tkscrollbar(w, command =
                               function(...) tcltk::tkyview(can, ...))
    tcltk::tkconfigure(can, yscrollcommand =
                       function(...) tcltk::tkset(yscr, ...))
    ## Pack it up
    tcltk::tkpack(buts, side="bottom", fill="x", pady="2m")
    tcltk::tkpack(can, side="left", fill="x")
    tcltk::tkpack(yscr, side="right", fill="y")
    tcltk::tkgrid(cp2eps, export, dump, dismiss, sticky="s")
    ## Box
    x0 <- usr2xy(c(xrange[1], yrange[1]))
    x1 <- usr2xy(c(xrange[2], yrange[2]))
    tcltk::tkcreate(can, "rectangle", x0[1], x0[2], x1[1], x1[2], outline = p$fg,
             width = p$lwd)
    ## Axes and ticks
    tl <-  -p$tcl * rpix     # -p$tcl * p$ps * p2p
    axoff <- p$mgp[3] * rpix
    tmp <- xpretty
    for (i in seq_along(tmp)) {
        x0 <- usr2xy(c(xpretty[1], yrange[1]))
        x1 <- usr2xy(c(xpretty[length(xpretty)], yrange[1]))
        tcltk::tkcreate(can, "line", x0[1], x0[2]+axoff, x1[1], x1[2]+axoff,
                 fill=p$fg)
        xx <- usr2xy(c(tmp[i], yrange[1]))
        tcltk::tkcreate(can, "line", xx[1], xx[2] + axoff, xx[1],
                        xx[2]+tl+axoff, fill=p$fg)
        tcltk::tkcreate(can, "text", xx[1], xx[2] + rpix * p$mgp[2], anchor="n",
                 text=as.character(tmp[i]), fill=p$col.axis, font=fnt.axis)
    }
    xx <- usr2xy(c(mean(xrange), yrange[1]))
    tcltk::tkcreate(can, "text", xx[1], xx[2] + rpix * p$mgp[1],
             text=colnames(sco)[1], fill=p$col.lab, anchor="n", font=fnt.lab)
    tmp <- ypretty
    for (i in seq_along(tmp)) {
        x0 <- usr2xy(c(xrange[1], tmp[1]))
        x1 <- usr2xy(c(xrange[1], tmp[length(tmp)]))
        tcltk::tkcreate(can, "line", x0[1]-axoff, x0[2], x1[1]-axoff, x1[2])
        yy <- usr2xy(c(xrange[1], tmp[i]))
        tcltk::tkcreate(can, "line", yy[1]-axoff, yy[2], yy[1]-tl-axoff, yy[2],
                 fill=p$fg )
        tcltk::tkcreate(can, "text", yy[1] - rpix * p$mgp[2] , yy[2], anchor="e",
                 text=as.character(tmp[i]), fill = p$col.axis, font=fnt.axis)
    }
    ## Points and labels

    ## The following 'inherits' works with ordipointlabel, but not
    ## with zooming
    if (inherits(x, "orditkplot")) {
        lsco <- scores(x, "labels")
        laboff <- rep(0, nrow(lsco))
        lsco <- lsco[rownames(sco),]
    } else {
        lsco <- sco
        laboff <- round(p2p * p$ps/2 + diam + 1)
    }
    pola <- tcltk::tclArray()        # points
    labtext <- tcltk::tclArray()     # text
    id <- tcltk::tclArray()          # index
    for (i in 1:nrow(sco)) {
        xy <- usr2xy(sco[i,])
        item <- Point(xy[1], xy[2], pch = pch[i], col = pcol[i],
                      fill = pbg[i], diam = diam[i])
        xy <- usr2xy(lsco[i,])
        fnt <- c(labfam[i], labsize[i], saneslant(labfnt[i]))
        lab <- tcltk::tkcreate(can, "text", xy[1], xy[2]-laboff[i], text=labs[i],
                        fill = tcol[i], font=fnt)
        tcltk::tkaddtag(can, "point", "withtag", item)
        tcltk::tkaddtag(can, "label", "withtag", lab)
        pola[[lab]] <- item
        labtext[[lab]] <- labs[i]
        id[[lab]] <- i
    }

##############################
### Mouse operations on canvas
##############################

    ## Plotting and Moving
    ## Mouse enters a label
    pEnter <- function() {
        tcltk::tkdelete(can, "box")
        hbox <- tcltk::tkcreate(can, "rectangle",
                                tcltk::tkbbox(can, "current"),
                         outline = "red", fill = "yellow")
        tcltk::tkaddtag(can, "box", "withtag", hbox)
        tcltk::tkitemraise(can, "current")
    }
    ## Mouse leaves a label
    pLeave <- function() {
        tcltk::tkdelete(can, "box")
    }
    ## Select label
    pDown <- function(x, y) {
        x <- as.numeric(x)
        y <- as.numeric(y)
        tcltk::tkdtag(can, "selected")
        tcltk::tkaddtag(can, "selected", "withtag", "current")
        tcltk::tkitemraise(can, "current")
        p <- as.numeric(tcltk::tkcoords(can,
                             pola[[tcltk::tkfind(can, "withtag", "current")]]))
        .pX <<- (p[1]+p[3])/2
        .pY <<- (p[2]+p[4])/2
        .lastX <<- x
        .lastY <<- y
    }
    ## Move label
    pMove <- function(x, y) {
        x <- as.numeric(x)
        y <- as.numeric(y)
        tcltk::tkmove(can, "selected", x - .lastX, y - .lastY)
        tcltk::tkdelete(can, "ptr")
        tcltk::tkdelete(can, "box")
        .lastX <<- x
        .lastY <<- y
        ## xadj,yadj: adjust for canvas scrolling
        xadj <- as.numeric(tcltk::tkcanvasx(can, 0))
        yadj <- as.numeric(tcltk::tkcanvasy(can, 0))
        hbox <- tcltk::tkcreate(can, "rectangle",
                                tcltk::tkbbox(can, "selected"),
                                outline = "red")
        tcltk::tkaddtag(can, "box", "withtag", hbox)
        conn <- tcltk::tkcreate(can, "line", .lastX + xadj, .lastY+yadj,
                         .pX, .pY, fill="red")
        tcltk::tkaddtag(can, "ptr", "withtag", conn)
    }
    ## Edit label
    pEdit <- function() {
        tcltk::tkdtag(can, "selected")
        tcltk::tkaddtag(can, "selected", "withtag", "current")
        tcltk::tkitemraise(can, "current")
        click <- tcltk::tkfind(can, "withtag", "current")
        txt <- tcltk::tclVar(labtext[[click]])
        i <- as.numeric(id[[click]])
        tt <- tcltk::tktoplevel()
        labEd <- tcltk::tkentry(tt, width=20, textvariable=txt)
        tcltk::tkgrid(tcltk::tklabel(tt, text = "Edit label"))
        tcltk::tkgrid(labEd, pady="5m", padx="5m")
        isDone <- function() {
            txt <- tcltk::tclvalue(txt)
            tcltk::tkitemconfigure(can, click, text = txt)
            rownames(sco)[i] <<- txt
            tcltk::tkdestroy(tt)
        }
        tcltk::tkbind(labEd, "<Return>", isDone)
    }
    ## Zooming: draw rectangle and take its user coordinates
    ## Rectangle: first corner
    pRect0 <- function(x, y) {
        x <- as.numeric(x)
        y <- as.numeric(y)
        ## yadj here and below adjusts for canvas scrolling
        yadj <- as.numeric(tcltk::tkcanvasy(can, 0))
        .pX <<- x
        .pY <<- y + yadj
    }
    ## Grow rectangle
    pRect <- function(x, y) {
        x <- as.numeric(x)
        y <- as.numeric(y)
        tcltk::tkdelete(can, "box")
        yadj <- as.numeric(tcltk::tkcanvasy(can, 0))
        .lastX <<- x
        .lastY <<- y + yadj
        rect <- tcltk::tkcreate(can, "rectangle", .pX, .pY, .lastX, .lastY,
                         outline="blue")
        tcltk::tkaddtag(can, "box", "withtag", rect)
    }
    ## Redraw ordiktplot with new xlim and ylim
    pZoom <- function() {
        nxlim <- sort(c(x2usr(.pX), x2usr(.lastX)))
        nylim <- sort(c(y2usr(.pY), y2usr(.lastY)))
        xy <- ordDump()
        ## Move labels closer to points in zoom
        ## FIXME: Doesn't do a perfect job
        mul <- abs(diff(nxlim)/diff(xlim))
        xy$labels <- xy$points + (xy$labels - xy$points)*mul
        xy$args$xlim <- nxlim
        xy$args$ylim <- nylim
        orditkplot(xy)
    }
    ## Dummy location of the mouse
    .lastX <- 0
    .lastY <- 0
    .pX <- 0
    .pY <- 0
    ## Mouse bindings:
    ## Moving a label
    tcltk::tkitembind(can, "label", "<Any-Enter>", pEnter)
    tcltk::tkitembind(can, "label", "<Any-Leave>", pLeave)
    tcltk::tkitembind(can, "label", "<1>", pDown)
    tcltk::tkitembind(can, "label", "<ButtonRelease-1>",
               function() {tcltk::tkdtag(can, "selected")
                           tcltk::tkdelete(can, "ptr")})
    tcltk::tkitembind(can, "label", "<B1-Motion>", pMove)
    ## Edit labels
    tcltk::tkitembind(can, "label", "<Double-Button-1>", pEdit)
    ## Zoom (with one-button mouse)
    tcltk::tkbind(can, "<Shift-Button-1>", pRect0)
    tcltk::tkbind(can, "<Shift-B1-Motion>", pRect)
    tcltk::tkbind(can, "<Shift-ButtonRelease>", pZoom)
    ## Zoom (with right button)
    tcltk::tkbind(can, "<Button-3>", pRect0)
    tcltk::tkbind(can, "<B3-Motion>", pRect)
    tcltk::tkbind(can, "<ButtonRelease-3>", pZoom)
}

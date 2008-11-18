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
    require(tcltk) || stop("requires package tcltk")

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
    p2p <- as.numeric(tclvalue(tcl("tk", "scaling"))) # Pixel per point
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
               "3" = {tkcreate(can, "line",
                               x, y+SQ*diam, x, y-SQ*diam, fill=col)
                      tkcreate(can, "line",
                               x+SQ*diam, y, x-SQ*diam, y, fill=col)},
               "4" = {tkcreate(can, "line",
                               x-diam, y-diam, x+diam, y+diam, fill=col)
                      tkcreate(can, "line",
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
               "14" = {tkcreate(can, "line", x-diam, y-diam, x, y+diam,
                                fill = col)
                       tkcreate(can, "line", x+diam, y-diam, x, y+diam,
                                fill = col)
                       Point(x, y, 0, col, fill, diam)},
               "15" = Point(x, y, 22, col = col, fill = col, diam),
               "16" = Point(x, y, 21, col = col, fill = col, diam),
               "17" = Point(x, y, 24, col = col, fill = col, diam),
               "18" = Point(x, y, 23, col = col, fill = col, diam/SQ),
               "19" = Point(x, y, 21, col = col, fill = col, diam),
               "20" = Point(x, y, 21, col = col, fill = col, diam/2),
               "21" = tkcreate(can, "oval", x-diam, y-diam,
               x+diam, y+diam, outline = col, fill = fill),
               "22" = tkcreate(can, "rectangle", x-diam, y-diam,
               x+diam, y+diam, outline = col, fill = fill),
               "23" = tkcreate(can, "polygon", x, y+SQ*diam,
               x+SQ*diam, y, x, y-SQ*diam, x-SQ*diam, y,
               outline = col, fill = fill),
               "24" = tkcreate(can, "polygon", x, y-SQ*diam,
               x+sqrt(6)/2*diam, y+SQ/2*diam, x-sqrt(6)/2*diam, y+SQ/2*diam,
               outline = col, fill = fill),
               "25" = tkcreate(can, "polygon", x, y+SQ*diam,
               x+sqrt(6)/2*diam, y-SQ/2*diam, x-sqrt(6)/2*diam, y-SQ/2*diam,
               outline = col, fill = fill),
               "o" = Point(x, y, 1, col, fill, diam),
               ## default: text with dummy location of the label
               {tkcreate(can, "text",
                        x, y, text = as.character(pch), fill = col)
                Point(x, y, 21, col="", fill="", diam)}
                )
    }

############################
### Initialize Tcl/Tk Window
############################

    ## toplevel
    w <- tktoplevel()
    tktitle(w) <- deparse(match.call())
    ## Max dim of windows (depends on screen)
    YSCR <- as.numeric(tkwinfo("screenheight", w)) - 150
    XSCR <- as.numeric(tkwinfo("screenwidth", w))

################################
### Buttons and button functions
################################

    ## Buttons
    buts <- tkframe(w)
    ## Copy current canvas to EPS using the standard Tcl/Tk utility
    cp2eps <- tkbutton(buts, text="Copy to EPS", 
                       command=function() tkpostscript(can, x=0, y=0,
                       height=height, width=width, 
                       file=tkgetSaveFile(filetypes="{{EPS file} {.eps}}")))
    dismiss <- tkbutton(buts, text="Dismiss", command=function() tkdestroy(w))
    ## Dump current plot to an "orditkplot" object (internally)
    ordDump <- function() {
        xy <- matrix(0, nrow=nrow(sco), ncol=2)
        rownames(xy) <- rownames(sco)
        colnames(xy) <- colnames(sco)
        for(nm in names(pola)) {
            xy[as.numeric(tclvalue(id[[nm]])),] <- xy2usr(nm)
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
        dumpVar <- tclVar("")
        tt <- tktoplevel()
        tktitle(tt) <- "R Dump"
        entryDump <- tkentry(tt, width=20, textvariable=dumpVar)
        tkgrid(tklabel(tt, text="Enter name for an R object"))
        tkgrid(entryDump, pady="5m")
        isDone <- function() {
            dumpName <- tclvalue(dumpVar)
            if (exists(dumpName, envir=.GlobalEnv)) {
                ok <- tkmessageBox(message=paste(sQuote(dumpName),
                                   "exists.\nOK to overwrite?"),
                                   icon="warning", type="okcancel",
                                   default="ok")
                if(tclvalue(ok) == "ok") {
                    assign(dumpName, xy, envir=.GlobalEnv)
                    tkdestroy(tt)
                }
            }
            else {
                assign(dumpName, xy, envir=.GlobalEnv)
                tkdestroy(tt)
            }
        }
        tkbind(entryDump, "<Return>", isDone)
        tkfocus(tt)
    }
    dump <- tkbutton(buts, text="Dump to R", command=pDump)
    ## Button to write current "orditkplot" object to a graphical device
    devDump <- function() {
        xy <- ordDump()
        ftypes <- c("eps" = "{EPS File} {.eps}",
                    "pdf" = "{PDF File} {.pdf}",
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
        ## Should work also in R < 2.8.0 with no capabilities("tiff")
        if (!isTRUE(unname(capabilities("tiff"))))
            falt["tiff"] <- FALSE
        ## bmp lives only in Windows
        if (.Platform$OS.type != "windows")
            falt["bmp"] <- FALSE
        ftypes <- ftypes[falt]
        fname <- tkgetSaveFile(filetypes=ftypes)
        if(tclvalue(fname) == "")
            return(NULL)
        fname <- tclvalue(fname)
        ftype <- unlist(strsplit(fname, "\\."))
        ftype <- ftype[length(ftype)]
        if (ftype == "jpeg")
            ftype <- "jpg"
        if (ftype == "tif")
            ftype <- "tiff"
        mess <- "is not a supported type: file not produced. Supported types are"
        if (!(ftype %in% names(ftypes))) {
            tkmessageBox(message=paste(sQuote(ftype), mess, paste(names(ftypes),
                         collapse=", ")), icon="warning")
            return(NULL)
        }
        pixdim <- round(xy$dim*PPI*p2p)
        switch(ftype,
               eps = postscript(file=fname, width=xy$dim[1], height=xy$dim[2],
               paper="special", horizontal = FALSE),
               pdf = pdf(file=fname, width=xy$dim[1], height=xy$dim[2]),
               png = png(file=fname, width=pixdim[1], height=pixdim[2]),
               jpg = jpeg(file=fname, width=pixdim[1], height=pixdim[2],
               quality = 100),
               tiff = tiff(file=fname, width=pixdim[1], height=pixdim[2]),
               bmp = bmp(file=fname, width=pixdim[1], height=pixdim[2]),
               fig = xfig(file=fname, width=xy$dim[1], height=xy$dim[2]))
        plot.orditkplot(xy)
        dev.off()
    }
    export <- tkbutton(buts, text="Export plot", command=devDump)

##########
### Canvas
##########
    
    ## Make canvas
    sco <- try(scores(x, display=display, choices = choices, ...),
               silent = TRUE)
    if (inherits(sco, "try-error")) {
        tkmessageBox(message=paste("No ordination scores were found in",
                     sQuote(deparse(substitute(x)))), icon="error")
        tkdestroy(w)
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
        xlim <- range(sco[,1])
    if (missing(ylim))
        ylim <- range(sco[,2])
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
        xy <- as.numeric(tkcoords(can, item))
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
    can <- tkcanvas(w, relief="sunken", width=width, height=min(height,YSCR),
                    scrollregion=c(0,0,width,height))
    if (p$bg != "")
        tkconfigure(can, bg=p$bg)
    yscr <- tkscrollbar(w, command = function(...) tkyview(can, ...))
    tkconfigure(can, yscrollcommand = function(...) tkset(yscr, ...))
    ## Pack it up
    tkpack(buts, side="bottom", fill="x", pady="2m")
    tkpack(can, side="left", fill="x")
    tkpack(yscr, side="right", fill="y")
    tkgrid(cp2eps, export, dump, dismiss, sticky="s")
    ## Box
    x0 <- usr2xy(c(xrange[1], yrange[1]))
    x1 <- usr2xy(c(xrange[2], yrange[2]))
    tkcreate(can, "rectangle", x0[1], x0[2], x1[1], x1[2], outline = p$fg,
             width = p$lwd)
    ## Axes and ticks
    tl <-  -p$tcl * rpix     # -p$tcl * p$ps * p2p
    axoff <- p$mgp[3] * rpix
    tmp <- xpretty
    for (i in 1:length(tmp)) {
        x0 <- usr2xy(c(xpretty[1], yrange[1]))
        x1 <- usr2xy(c(xpretty[length(xpretty)], yrange[1]))
        tkcreate(can, "line", x0[1], x0[2]+axoff, x1[1], x1[2]+axoff,
                 fill=p$fg)
        xx <- usr2xy(c(tmp[i], yrange[1]))
        tkcreate(can, "line", xx[1], xx[2] + axoff, xx[1], xx[2]+tl+axoff,
                 fill=p$fg)
        tkcreate(can, "text", xx[1], xx[2] + rpix * p$mgp[2], anchor="n",
                 text=as.character(tmp[i]), fill=p$col.axis, font=fnt.axis)
    }
    xx <- usr2xy(c(mean(xrange), yrange[1]))
    tkcreate(can, "text", xx[1], xx[2] + rpix * p$mgp[1],
             text=colnames(sco)[1], fill=p$col.lab, anchor="n", font=fnt.lab)
    tmp <- ypretty
    for (i in 1:length(tmp)) {
        x0 <- usr2xy(c(xrange[1], tmp[1]))
        x1 <- usr2xy(c(xrange[1], tmp[length(tmp)]))
        tkcreate(can, "line", x0[1]-axoff, x0[2], x1[1]-axoff, x1[2])
        yy <- usr2xy(c(xrange[1], tmp[i]))
        tkcreate(can, "line", yy[1]-axoff, yy[2], yy[1]-tl-axoff, yy[2],
                 fill=p$fg )
        tkcreate(can, "text", yy[1] - rpix * p$mgp[2] , yy[2], anchor="e",
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
    pola <- tclArray()        # points
    labtext <- tclArray()     # text
    id <- tclArray()          # index
    for (i in 1:nrow(sco)) {
        xy <- usr2xy(sco[i,])
        item <- Point(xy[1], xy[2], pch = pch[i], col = pcol[i],
                      fill = pbg[i], diam = diam[i])
        xy <- usr2xy(lsco[i,])
        fnt <- c(labfam[i], labsize[i], saneslant(labfnt[i]))
        lab <- tkcreate(can, "text", xy[1], xy[2]-laboff[i], text=labs[i],
                        fill = tcol[i], font=fnt)
        tkaddtag(can, "point", "withtag", item)
        tkaddtag(can, "label", "withtag", lab)
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
        tkdelete(can, "box")
        hbox <- tkcreate(can, "rectangle", tkbbox(can, "current"),
                         outline = "red", fill = "yellow")
        tkaddtag(can, "box", "withtag", hbox)
        tkitemraise(can, "current")
    }
    ## Mouse leaves a label
    pLeave <- function() {
        tkdelete(can, "box")
    }
    ## Select label
    pDown <- function(x, y) {
        x <- as.numeric(x)
        y <- as.numeric(y)
        tkdtag(can, "selected")
        tkaddtag(can, "selected", "withtag", "current")
        tkitemraise(can, "current")
        p <- as.numeric(tkcoords(can,
                                 pola[[tkfind(can, "withtag", "current")]]))
        .pX <<- (p[1]+p[3])/2
        .pY <<- (p[2]+p[4])/2
        .lastX <<- x
        .lastY <<- y
    }
    ## Move label
    pMove <- function(x, y) {
        x <- as.numeric(x)
        y <- as.numeric(y)
        tkmove(can, "selected", x - .lastX, y - .lastY)
        tkdelete(can, "ptr")
        tkdelete(can, "box")
        .lastX <<- x
        .lastY <<- y
        ## xadj,yadj: adjust for canvas scrolling
        xadj <- as.numeric(tkcanvasx(can, 0))
        yadj <- as.numeric(tkcanvasy(can, 0))
        hbox <- tkcreate(can, "rectangle", tkbbox(can, "selected"),
                         outline = "red")
        tkaddtag(can, "box", "withtag", hbox)
        conn <- tkcreate(can, "line", .lastX + xadj, .lastY+yadj,
                         .pX, .pY, fill="red")
        tkaddtag(can, "ptr", "withtag", conn)
    }
    ## Edit label
    pEdit <- function() {
        tkdtag(can, "selected")
        tkaddtag(can, "selected", "withtag", "current")
        tkitemraise(can, "current")
        click <- tkfind(can, "withtag", "current")
        txt <- tclVar(labtext[[click]])
        i <- as.numeric(id[[click]])
        tt <- tktoplevel()
        labEd <- tkentry(tt, width=20, textvariable=txt)
        tkgrid(tklabel(tt, text = "Edit label"))
        tkgrid(labEd, pady="5m", padx="5m")
        isDone <- function() {
            txt <- tclvalue(txt)
            tkitemconfigure(can, click, text = txt)
            rownames(sco)[i] <<- txt
            tkdestroy(tt)
        }
        tkbind(labEd, "<Return>", isDone)
    }
    ## Zooming: draw rectangle and take its user coordinates
    ## Rectangle: first corner
    pRect0 <- function(x, y) {
        x <- as.numeric(x)
        y <- as.numeric(y)
        ## yadj here and below adjusts for canvas scrolling
        yadj <- as.numeric(tkcanvasy(can, 0))
        .pX <<- x
        .pY <<- y + yadj
    }
    ## Grow rectangle
    pRect <- function(x, y) {
        x <- as.numeric(x)
        y <- as.numeric(y)
        tkdelete(can, "box")
        yadj <- as.numeric(tkcanvasy(can, 0))
        .lastX <<- x
        .lastY <<- y + yadj
        rect <- tkcreate(can, "rectangle", .pX, .pY, .lastX, .lastY,
                         outline="blue")
        tkaddtag(can, "box", "withtag", rect)
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
    tkitembind(can, "label", "<Any-Enter>", pEnter)
    tkitembind(can, "label", "<Any-Leave>", pLeave)
    tkitembind(can, "label", "<1>", pDown)
    tkitembind(can, "label", "<ButtonRelease-1>",
               function() {tkdtag(can, "selected"); tkdelete(can, "ptr")})
    tkitembind(can, "label", "<B1-Motion>", pMove)
    ## Edit labels
    tkitembind(can, "label", "<Double-Button-1>", pEdit) 
    ## Zoom (with one-button mouse)
    tkbind(can, "<Shift-Button-1>", pRect0)
    tkbind(can, "<Shift-B1-Motion>", pRect)
    tkbind(can, "<Shift-ButtonRelease>", pZoom)
    ## Zoom (with right button)
    tkbind(can, "<Button-3>", pRect0)
    tkbind(can, "<B3-Motion>", pRect)
    tkbind(can, "<ButtonRelease-3>", pZoom)
}

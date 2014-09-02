`ordisurf` <-
    function(...) UseMethod("ordisurf")

`ordisurf.formula` <-
    function(formula, data, ...)
{
    if (missing(data))
        data <- parent.frame()
    x <- formula[[2]]
    x <- eval.parent(x)
    formula[[2]] <- NULL
    y <- drop(as.matrix(model.frame(formula, data, na.action = na.pass)))
    if (NCOL(y) > 1)
        stop(gettextf("only one fitted variable allowed in the formula"))
    ordisurf(x, y, ...)
}

`ordisurf.default` <-
    function (x, y, choices = c(1, 2), knots = 10, family = "gaussian",
              col = "red", isotropic = TRUE, thinplate = TRUE, bs = "tp",
              fx = FALSE, add = FALSE, display = "sites", w = weights(x),
              main, nlevels = 10, levels, npoints = 31, labcex = 0.6,
              bubble = FALSE, cex = 1, select = TRUE, method = "REML",
              gamma = 1, plot = TRUE, lwd.cl = par("lwd"), ...)
{
    weights.default <- function(object, ...) NULL
    if(!missing(thinplate)) {
        warning("use of 'thinplate' is deprecated and will soon be removed;\nuse 'isotropic' instead")
        isotropic <- thinplate
    }
    ## GRID no user-definable - why 31?
    GRID <- npoints
    w <- eval(w)
    if (!is.null(w) && length(w) == 1)
        w <- NULL
    X <- scores(x, choices = choices, display = display, ...)
    ## The original name of 'y' may be lost in handling NA: save for
    ## plots
    yname <- deparse(substitute(y))
    kk <- complete.cases(X) & !is.na(y)
    if (!all(kk)) {
        X <- X[kk, , drop = FALSE]
        y <- y[kk]
        w <- w[kk]
    }
    x1 <- X[, 1]
    x2 <- X[, 2]
    ## handle fx - allow vector of length up to two
    if(!(missfx <- missing(fx)) && missing(knots))
        warning("requested fixed d.f. splines but without specifying 'knots':\nswitching to 'fx = FALSE'")
    if (length(fx) > 2L)
        warning("length of 'fx' supplied exceeds '2': using the first two")
    ## expand fx robustly, no matter what length supplied
    fx <- rep(fx, length.out = 2)
    ## can't have `fx = TRUE` and `select = TRUE`
    if(!missfx) { ## fx set by user
        if((miss.select <- missing(select)) && any(fx)) {
            warning("'fx = TRUE' requested; using 'select = FALSE'")
            select <- FALSE
        } else if(!miss.select && isTRUE(select)){
            stop("fixed d.f. splines ('fx = TRUE') incompatible with 'select = TRUE'")
        }
    }
    ## handle knots - allow vector of length up to two
    if (length(knots) > 2L)
        warning("length of 'knots' supplied exceeds '2': using the first two")
    ## expand knots robustly, no matter what length supplied
    knots <- rep(knots, length.out = 2)
    ## handle the bs - we only allow some of the possible options
    if (length(bs) > 2L)
        warning("number of basis types supplied exceeds '2': only using the first two")
    bs <- rep(bs, length.out = 2)
    ## check allowed types
    BS <- c("tp","ts","cr","cs","ds","ps","ad")
    want <- match(bs, BS)
    user.bs <- bs ## store supplied (well expanded supplied ones)
    bs <- BS[want]
    if (any(wrong <- is.na(bs))) {
        stop(paste("Supplied basis type of",
                   paste(sQuote(unique(user.bs[wrong])), collapse = ", "),
                   "not supported."))
    }
    ## can't use "cr", "cs", "ps" in 2-d smoother with s()
    if(isTRUE(isotropic) && any(bs %in% c("cr", "cs", "ps"))) {
        stop("bases \"cr\", \"cs\", and \"ps\" not allowed in isotropic smooths")
    }
    ## Build formula
    if (knots[1] <= 0) {
        f <- formula(y ~ x1 + x2)
    } else if (knots[1] == 1) { ## why do we treat this differently?
        f <- formula(y ~ poly(x1, 1) + poly(x2, 1))
    } else if (knots[1] == 2) {
        f <- formula(y ~ poly(x1, 2) + poly(x2, 2) + poly(x1, 1):poly(x2, 1))
    } else if (isotropic) {
        f <- formula(paste0("y ~ s(x1, x2, k = ", knots[1],
                            ", bs = \"", bs[1], "\", fx = ", fx[1],")"))
    } else {
        if (any(bs %in% c("ad"))) {
            ## only "ad" for now, but "fs" should also not be allowed
            f <- formula(paste0("y ~ s(x1, k = ", knots[1],
                                ", bs = \"", bs[1],
                                "\", fx = ", fx[1], ") + s(x2, k = ",
                                knots[2], ", bs = \"", bs[2],
                                "\", fx = ", fx[2], ")"))
        } else {
            f <- formula(paste0("y ~ te(x1, x2, k = c(",
                                paste0(knots, collapse = ", "),
                                "), bs = c(",
                                paste0("\"", bs, "\"", collapse = ", "),
                                "), fx = c(",
                                paste0(fx, collapse = ", "),
                                "))"))
        }
    }
    ## fit model
    mod <- gam(f, family = family, weights = w, select = select,
               method = method, gamma = gamma)
    xn1 <- seq(min(x1), max(x1), len=GRID)
    xn2 <- seq(min(x2), max(x2), len=GRID)
    newd <- expand.grid(x1 = xn1, x2 = xn2)
    fit <- predict(mod, type = "response", newdata=as.data.frame(newd))
    poly <- chull(cbind(x1,x2))
    ## Move out points of the convex hull to have contour for all data
    ## points
    xhull1 <- x1[poly] + sign(x1[poly] - mean(x1[poly])) *
        diff(range(x1))/(GRID - 1)
    xhull2 <- x2[poly] + sign(x2[poly] - mean(x2[poly])) *
        diff(range(x2))/(GRID - 1)
    npol <- length(poly)
    np <- nrow(newd)
    inpoly <- numeric(np)
    inpoly <- .C("pnpoly", as.integer(npol), as.double(xhull1),
                 as.double(xhull2), as.integer(np), as.double(newd[,1]),
                 as.double(newd[,2]), inpoly = as.integer(inpoly),
                 PACKAGE="vegan")$inpoly
    is.na(fit) <- inpoly == 0
    if(plot) {
        if (!add) {
            if (bubble) {
                if (is.numeric(bubble))
                    cex <- bubble
                cex <- (y -  min(y))/diff(range(y)) * (cex-0.4) + 0.4
            }
            plot(X, asp = 1, cex = cex, ...)
        }
        if (!missing(main) || (missing(main) && !add)) {
            if (missing(main))
                main <- yname
            title(main = main)
        }
        if (missing(levels))
            levels <- pretty(range(fit, finite = TRUE), nlevels)
        ## Only plot surface is select is FALSE or (TRUE and EDF is diff from 0)
        if(!select ||
           (select && !isTRUE(all.equal(as.numeric(summary(mod)$edf), 0))))
            contour(xn1, xn2, matrix(fit, nrow=GRID), col = col, add = TRUE,
                    levels = levels, labcex = labcex,
                    drawlabels = !is.null(labcex) && labcex > 0,
                    lwd = lwd.cl)
    }
    mod$grid <- list(x = xn1, y = xn2, z = matrix(fit, nrow = GRID))
    class(mod) <- c("ordisurf", class(mod))
    mod
}

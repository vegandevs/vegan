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
    y <- as.matrix(model.frame(formula, data, na.action = na.pass))
    ordisurf(x, y, ...)
}

`ordisurf.default` <-
    function (x, y, choices = c(1, 2), what = c("contour", "surface"),
              knots = 10, family = "gaussian", col = "red", isotropic = TRUE,
              bs = "tp", fx = FALSE, add = FALSE, display = "sites", w, main,
              nlevels = 10, levels, npoints = 51, labcex = 0.6, bubble = FALSE,
              cex = 1, select = TRUE, method = "REML", gamma = 1, plot = TRUE,
              lwd.cl = par("lwd"), ...)
{
    if (NCOL(y) > 1)
        stop(gettextf("only one fitted variable allowed in the formula"))
    what <- match.arg(what)
    ## GRID no user-definable - why 31?
    GRID <- npoints
    if (missing(w))
        w <- if(is.atomic(x)) attr(x, "weights")
             else weights(x, display = display)
    if (!is.null(w) && length(w) == 1)
        w <- NULL
    X <- scores(x, choices = choices, display = display, ...)
    ## The original names of 'y' may be lost in handling NA: save for
    ## plots
    if (missing(main))
        main <- colnames(y) %||% deparse(substitute(y))
    axlabs <- colnames(X)
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
        warning("number of basis types supplied exceeds '2': using the first two")
    bs <- rep(bs, length.out = 2)
    ## check allowed types
    BS <- c("tp","ts","cr","cs","ds","ps","ad")
    want <- match(bs, BS)
    user.bs <- bs ## store supplied (well expanded supplied ones)
    bs <- BS[want]
    if (any(wrong <- is.na(bs))) {
        stop(gettextf("supplied basis type of '%s' not supported",
                   paste(unique(user.bs[wrong]), collapse = ", ")))
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
    ## grid for prediction
    xpand <- diff(range(x1, x2)) / 25 # 4% expansion like in axes
    xn1 <- seq(min(x1) - xpand, max(x1) + xpand, len=GRID)
    xn2 <- seq(min(x2) - xpand, max(x2) + xpand, len=GRID)
    newd <- expand.grid(x1 = xn1, x2 = xn2)
    fit <- predict(mod, type = "response", newdata=as.data.frame(newd))
    poly <- chull(cbind(x1,x2))
    xhull1 <- x1[poly] + sign(x1[poly] - mean(x1[poly])) * xpand
    xhull2 <- x2[poly] + sign(x2[poly] - mean(x2[poly])) * xpand
    npol <- length(poly)
    np <- nrow(newd)
    inpoly <- numeric(np)
    inpoly <- .C(pnpoly, as.integer(npol), as.double(xhull1),
                 as.double(xhull2), as.integer(np), as.double(newd[,1]),
                 as.double(newd[,2]), inpoly = as.integer(inpoly),
                 PACKAGE = "vegan")$inpoly
    is.na(fit) <- inpoly == 0
    ## poly() models do not save args
    if (is.null(mod$model$x1) || is.null(mod$model$x2)) {
        mod$model$x1 <- x1
        mod$model$x2 <- x2
    }
    mod$grid <- list(x = xn1, y = xn2, z = matrix(fit, nrow = GRID))
    class(mod) <- c("ordisurf", class(mod))
    if (plot) {
        ## Do not draw contours or surface if select is TRUE and gam
        ## fit has EDF nearly 0 (fixed knots 0, 1 or 2 have no EDF,
        ## hence isTRUE)
        if (!select || !isTRUE(summary(mod)$edf < 1e-4))
            plot(mod, add = add, what = what, col = col, bubble = bubble,
                 cex = cex, nlevels = nlevels, levels = levels, labcex = labcex,
                 lwd.cl = lwd.cl, xlab = axlabs[1], ylab = axlabs[2],
                 main = main, ...)
    }
    mod
}

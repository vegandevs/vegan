`decostand` <-
    function (x, method, MARGIN, range.global, logbase = 2, na.rm = FALSE, ...)
{
    wasDataFrame <- is.data.frame(x)
    x <- as.matrix(x)
    METHODS <- c("total", "max", "frequency", "normalize", "range", "rank",
                 "rrank", "standardize", "pa", "chi.square", "hellinger",
                 "log", "clr", "rclr", "alr")
    method <- match.arg(method, METHODS)
    if (any(x < 0, na.rm = TRUE)) {
        k <- min(x, na.rm = TRUE)
        if (method %in% c("total", "frequency", "pa", "chi.square", "rank",
                          "rrank", "rclr")) {
            warning("input data contains negative entries: result may be non-sense")
        }
    }
    else k <- .Machine$double.eps
    attr <- NULL
    switch(method, total = {
        if (missing(MARGIN))
            MARGIN <- 1
        tmp <- pmax.int(k, apply(x, MARGIN, sum, na.rm = na.rm))
        x <- sweep(x, MARGIN, tmp, "/")
        attr <- list("total" = tmp, "margin" = MARGIN)
    }, max = {
        if (missing(MARGIN))
            MARGIN <- 2
        tmp <- pmax.int(k, apply(x, MARGIN, max, na.rm = na.rm))
        x <- sweep(x, MARGIN, tmp, "/")
        attr <- list("max" = tmp, "margin" = MARGIN)
    }, frequency = {
        if (missing(MARGIN))
            MARGIN <- 2
        tmp <- pmax.int(k, apply(x, MARGIN, sum, na.rm = na.rm))
        fre <- apply(x > 0, MARGIN, sum, na.rm = na.rm)
        tmp <- fre/tmp
        x <- sweep(x, MARGIN, tmp, "*")
        attr <- list("scale" = tmp, "margin" = MARGIN)
    }, normalize = {
        if (missing(MARGIN))
            MARGIN <- 1
        tmp <- apply(x^2, MARGIN, sum, na.rm = na.rm)
        tmp <- pmax.int(.Machine$double.eps, sqrt(tmp))
        x <- sweep(x, MARGIN, tmp, "/")
        attr <- list("norm" = tmp, "margin" = MARGIN)
    }, range = {
        if (missing(MARGIN))
            MARGIN <- 2
        if (missing(range.global))
            xtmp <- x
        else {
            if (dim(range.global)[MARGIN] != dim(x)[MARGIN])
                stop("range matrix does not match data matrix")
            xtmp <- as.matrix(range.global)
        }
        tmp <- apply(xtmp, MARGIN, min, na.rm = na.rm)
        ran <- apply(xtmp, MARGIN, max, na.rm = na.rm)
        ran <- ran - tmp
        ran <- pmax.int(k, ran, na.rm = na.rm)
        x <- sweep(x, MARGIN, tmp, "-")
        x <- sweep(x, MARGIN, ran, "/")
        attr <- list("min" = tmp, "range" = ran, "margin" = MARGIN)
    }, rank = {
        wasNA <- is.na(x)
        if (any(wasNA) && !na.rm)
            stop("missing values are not allowed with 'na.rm = FALSE'")
        if (missing(MARGIN)) MARGIN <- 1
        x[x==0] <- NA
        x <- apply(x, MARGIN, rank, na.last = "keep")
        if (MARGIN == 1) # gives transposed x
            x <- t(x)
        x[is.na(x)] <- 0
        if(any(wasNA))
            x[wasNA] <- NA
        attr <- list("margin" = MARGIN)
    }, rrank = {
        if (missing(MARGIN)) MARGIN <- 1
        x <- decostand(x, "rank", MARGIN = MARGIN, na.rm = na.rm)
        if (na.rm && any(wasNA <- is.na(x)))
            x[wasNA] <- 0
        x <- sweep(x, MARGIN, specnumber(x, MARGIN = MARGIN), "/")
        if (any(wasNA))
            x[wasNA] <- NA
        attr <- list("margin" = MARGIN)
    }, standardize = {
        if (!missing(MARGIN) && MARGIN == 1)
            x <- t(scale(t(x)))
        else {
            x <- scale(x)
            MARGIN <- 2
        }
        attr <- list("center" = attr(x, "scaled:center"),
                     "scale" = attr(x, "scaled:scale"),
                     "margin" = MARGIN)
    }, pa = {
        x <- ifelse(x > 0, 1, 0)
    }, chi.square = {
        if (missing(MARGIN))
            MARGIN <- 1
        ## MARGIN 2 transposes the result!
        if (MARGIN == 2)
            x <- t(x)
        rs <- pmax.int(k, rowSums(x, na.rm = na.rm))
        cs <- pmax.int(k, colSums(x, na.rm = na.rm))
        tot <- sum(x, na.rm = na.rm)
        x <- sqrt(tot) * x/outer(rs, sqrt(cs))
        attr <- list("tot" = tot, "rsum" = rs, "csum" = cs, margin = MARGIN)
    }, hellinger = {
        x <- sqrt(decostand(x, "total", MARGIN = MARGIN, na.rm = na.rm))
        attr <- attr(x, "parameters")
    }, log = {### Marti Anderson logs, after Etienne Laliberte
        if (!isTRUE(all.equal(as.integer(x), as.vector(x)))) {
            minpos <- min(x[x > 0], na.rm = TRUE)
            x <- x / minpos
            warning("non-integer data: divided by smallest positive value",
                    call. = FALSE)
        } else {
            minpos <- 1
        }
        x[x > 0 & !is.na(x)] <- log(x[x > 0 & !is.na(x)], base = logbase) + 1
        attr <- list("logbase" = logbase, minpos = minpos)
    }, alr = {
        if (missing(MARGIN))
	    MARGIN <- 1
        if (MARGIN == 1)
            x <- .calc_alr(x, na.rm=na.rm, ...)
	else x <- t(.calc_alr(t(x), na.rm=na.rm, ...))
        attr <- attr(x, "parameters")
        attr$margin <- MARGIN
    }, clr = {
        if (missing(MARGIN))
	    MARGIN <- 1
        if (MARGIN == 1)
            x <- .calc_clr(x, na.rm=na.rm, ...)
	else x <- t(.calc_clr(t(x), na.rm=na.rm, ...))
        attr <- attr(x, "parameters")
        attr$margin <- MARGIN
    }, rclr = {
        if (missing(MARGIN))
	    MARGIN <- 1
        if (MARGIN == 1)
            x <- .calc_rclr(x, ...)
	else x <- t(.calc_rclr(t(x), ...))
        attr <- attr(x, "parameters")
        attr$margin <- MARGIN
    })
    if (any(is.nan(x)))
        warning("result contains NaN, perhaps due to impossible mathematical
                 operation\n")
    if (wasDataFrame)
        x <- as.data.frame(x)
    attr(x, "parameters") <- attr
    attr(x, "decostand") <- method
    x
}

## Modified from the original version in mia R package
.calc_clr <-
    function(x, na.rm, pseudocount=0, ...)
{

    # Add pseudocount
    x <- x + pseudocount
    # Error with negative values (note: at this step we always use na.rm=TRUE)
    if (any(x <= 0, na.rm = TRUE)) {
        stop("'clr' cannot be used with non-positive data: use pseudocount > ",
             -min(x, na.rm = TRUE) + pseudocount, call. = FALSE)
    }

    # In every sample, calculate the log of individual entries.
    # Then calculate
    # the sample-specific mean value and subtract every entries'
    # value with that.
    clog <- log(x)

    # Calculate sample-wise log means (note: here we always set na.rm=TRUE !)
    means <- rowMeans(clog, na.rm = TRUE)

    # CLR transformation
    clog <- clog - means

    # Replace missing values with 0
    if (na.rm && any(is.na(clog))) {
        cat("Replacing missing values with zero for clr. You can disable this with na.rm=FALSE.")    
        clog[is.na(clog)] <- 0
    } 
    
    attr(clog, "parameters") <- list("means" = means,
                                     "pseudocount" = pseudocount)
    clog
}

# Modified from the original version in mia R package
.calc_rclr <-
    function(x, ROPT=3, NITER=5, TOL=1e-5, verbose=FALSE, impute=TRUE, ...)
{
    # Error with negative values
    if (any(x < 0)) {
        stop("'rclr' cannot be used with negative data", call. = FALSE)
    }

   # Log transform
   clog <- log(x)
   
   # Convert zeros to NAs in rclr
   clog[is.infinite(clog)] <- NA

   # Calculate log of geometric mean for every sample, ignoring the NAs
   # Always na.rm=TRUE at this step!   
   means <- rowMeans(clog, na.rm = TRUE)

   # TODO
   #if (any(is.na(means))) {
   #  stop("some samples do not contain (any) observations and thus rclr cannot be calculated")
   #}

   ## Divide (or in log-space, reduce) all values by their sample-wide
   ## geometric means
   xx <- clog - means
   attr(xx, "parameters") <- list("means" = means)

   # Impute NAs if impute=TRUE
   # Otherwise return the transformation with NAs
   if (impute && any(is.na(xx))){
   
     opt_res <- OptSpace(xx, ROPT=ROPT, NITER=NITER, TOL=TOL, verbose=verbose)
     
     # recenter the data (the means of rclr can get thrown off since we work on only missing)
     M_I <- opt_res$X %*% opt_res$S %*% t(opt_res$Y)

     # Center cols to 0
     M_I <- as.matrix(scale(M_I, center=TRUE, scale=FALSE))

     # Center rows to 0
     M_I <- as.matrix(t(scale(t(M_I), center=TRUE, scale=FALSE)))

     # Imputed matrix
     xx <- M_I
   }
  
   return(xx)
   
}

.calc_alr <-
    function (x, na.rm, pseudocount = 0, reference = 1, ...)
{
    # Add pseudocount
    x <- x + pseudocount
    # If there is negative values, gives an error.
    # Always na.rm=TRUE at this step
    if (any(x <= 0, na.rm = TRUE)) {
        stop("'alr' cannot be used with non-positive data: use pseudocount > ",
             -min(x, na.rm = na.rm) + pseudocount, call. = FALSE)
    }
    ## name must be changed to numeric index for [-reference,] to work
    if (is.character(reference)) {
        reference <- which(reference == colnames(x))
        if (!length(reference)) # found it?
            stop("'reference' name was not found in data", call. = FALSE)
    }
    if (reference > ncol(x) || reference < 1)
        stop("'reference' should be a name or index 1 to ",
             ncol(x), call. = FALSE)
    clog <- log(x)
    refvector <- clog[, reference]
    clog <- clog[, -reference] - refvector

    # Replace missing values with 0
    if (na.rm && any(is.na(clog))) {
        cat("Replacing missing values with zero for alr. You can disable this with na.rm=FALSE.")
        clog[is.na(clog)] <- 0
    } 

    attr(clog, "parameters") <- list("reference" = refvector,
                                     "index" = reference,
                                     "pseudocount" = pseudocount)
    clog
}

`decobackstand` <-
    function(x, zap = TRUE)
{
    method <- attr(x, "decostand")
    if (is.null(method))
        stop("function can be used only with 'decostand' standardized data")
    para <- attr(x, "parameters")
    if(is.null(para)) # for old results & "pa"
        stop("object has no information to backtransform data")
    x <- switch(method,
                "total" = sweep(x, para$margin, para$total, "*"),
                "max" = sweep(x, para$margin, para$max, "*"),
                "frequency" = sweep(x, para$margin, para$scale, "/"),
                "normalize" = sweep(x, para$margin, para$norm, "*"),
                "range" = { x <- sweep(x, para$margin, para$range, "*")
                            sweep(x, para$margin, para$min, "+")},
                "standardize" = {x <- sweep(x, para$margin, para$scale, "*")
                                 sweep(x, para$margin, para$center, "+") },
                "hellinger" = sweep(x^2, para$margin, para$total, "*"),
                "chi.square" = { rc <- outer(para$rsum, sqrt(para$csum))
                                 x <- x * rc /sqrt(para$tot)
                                 if (para$margin == 1) x else t(x) },
                "log" = { x[x > 0 & !is.na(x)] <-
                              para$logbase^(x[x > 0 & !is.na(x)] - 1)
                              x * para$minpos},
                "clr" = exp(sweep(x, para$margin, para$means, "+")) -
                    para$pseudocount,
                "rclr" = { x[x == 0] <- -Inf # x==0 was set: should be safe
                           exp(sweep(x, para$margin, para$means, "+"))},
                "wisconsin" = { x <- sweep(x, 1, para$total, "*")
                                sweep(x, 2, para$max, "*") },
                stop("no backtransformation available for method ",
                     sQuote(method))
                )
    if (zap)
        x[abs(x) < sqrt(.Machine$double.eps)] <- 0
    x
}


OptSpace <- function(x, ROPT=3, NITER=5, TOL=1e-5, verbose=FALSE)
{

  ## Preprocessing : x     : partially revealed matrix
  if (is.data.frame(x)) {
    x <- as.matrix(x)
  }
  if (!is.matrix(x)){
    stop("* OptSpace : an input x should be a matrix")
  }
  if (any(is.infinite(x))){
    stop("* OptSpace : no infinite value in x is allowed")
  }
  if (!any(is.na(x))){
    stop("* OptSpace : there is no unobserved values as NA")
  }
  idxna <- (is.na(x))
  M_E <- array(0, c(nrow(x), ncol(x)))
  M_E[!idxna] <- x[!idxna]
  
  ## Preprocessing : size information
  n <- nrow(x)
  m <- ncol(x)
  
  ## Preprocessing : other sparse-related concepts
  nnZ.E <- sum(!idxna)
  E <- array(0, c(nrow(x), ncol(x))); E[!idxna] <- 1
  eps <- nnZ.E/sqrt(m*n)
  
  ## Preprocessing : ROPT  : implied rank
  if (is.na(ROPT)){
    if (verbose){
      cat("* OptSpace: Guessing an implicit rank.")
    }
    r <- min(max(round(.guess_rank(M_E, nnZ.E)), 2), m-1)
    if (verbose){
      cat(paste0('* OptSpace: Guessing an implicit rank: Estimated rank : ',r))
    }
  } else {
    r <- round(ROPT)
    if ((!is.numeric(r)) || (r<1) || (r>m) || (r>n)){
      stop("* OptSpace: ROPT should be an integer in [1,min(nrow(x),ncol(x))]")
    }
  }
  
  ## Preprocessing : NITER : maximum number of iterations
  if ((is.infinite(NITER))||(NITER<=1)||(!is.numeric(NITER))){
    stop("* OptSpace: invalid NITER number")
  }
  NITER <- round(NITER)
  
  m0 <- 10000
  rho <-  eps*n
  
  ## Main Computation
  rescal_param <- sqrt(nnZ.E*r/(norm(M_E,'f')^2))
  M_E <- M_E*rescal_param
  
  # 1. SVD
  if (verbose){
    cat("* OptSpace: Step 2: SVD ...")
  }
  svdEt <- svd(M_E)
  X0 <- svdEt$u[,1:r]
  X0 <- X0[, rev(seq_len(ncol(X0)))]
  S0 <- diag(rev(svdEt$d[seq_len(r)]))
  Y0 <- svdEt$v[, seq_len(r)]
  Y0 <- Y0[, rev(seq_len(ncol(Y0)))]
  
  # 3. Initial Guess
  if (verbose){
    cat("* OptSpace: Step 3: Initial Guess ...")
  }
  X0 <- X0*sqrt(n)
  Y0 <- Y0*sqrt(m)
  S0 <- S0/eps
  
  # 4. Gradient Descent
  if (verbose){
    cat("* OptSpace: Step 4: Gradient Descent ...")
  }
  X <- X0
  Y <- Y0
  S <- .aux_getoptS(X, Y, M_E, E)
  # initialize
  dist <- array(0, c(1, (NITER+1)))
  dist[1] <- norm((M_E - (X %*% S %*% t(Y)))*E, 'f') / sqrt(nnZ.E)
  for (i in seq_len(NITER)){
    # compute the gradient
    tmpgrad <- .aux_gradF_t(X, Y, S, M_E, E, m0, rho)
    W <- tmpgrad$W
    Z <- tmpgrad$Z
    # line search for the optimum jump length
    t <- .aux_getoptT(X, W, Y, Z, S, M_E, E, m0, rho)
    X <- X + t*W;
    Y <- Y + t*Z;
    S <- .aux_getoptS(X, Y, M_E, E)
    # compute the distortion
    dist[i+1] <- norm(((M_E - X %*% S %*% t(Y))*E),'f') / sqrt(nnZ.E)
    if (dist[i+1] < TOL){
      dist <- dist[1:(i+1)]
      break
    }
  }
  S <- S/rescal_param  
  # Return Results
  out <- list()
  
  # re-order Optspace may change order during iters
  index_order <- order(diag(S), decreasing = TRUE)
  X <- X[, index_order]
  Y <- Y[, index_order]
  S <- S[index_order, index_order]
  out$X <- X
  out$S <- S
  out$Y <- Y
  out$dist <- dist
  if (verbose){
    cat('* OptSpace: estimation finished.')
  }
  return(out)
}


# @keywords internal
.guess_rank <- function(X, nnz)
{
  maxiter <- 10000
  n <- nrow(X)
  m <- ncol(X)
  epsilon <- nnz/sqrt(m*n)
  svdX <- svd(X)
  S0 <- svdX$d
  
  nsval0 <- length(S0)
  S1 <- S0[seq_len(nsval0-1)]-S0[seq(2, nsval0)]  
  nsval1 <- length(S1)
  if (nsval1 > 10){
    S1_ <- S1/mean(S1[seq((nsval1-10), nsval1)])
  } else {
    S1_ <- S1/mean(S1[seq_len(nsval1)])
  }
  r1 <- 0
  lam <- 0.05
  
  itcounter <- 0
  while (r1<=0){
    itcounter <- itcounter+1
    cost <- array(0, c(1, length(S1_)))
    for (idx in seq_len(length(S1_))) {
      cost[idx] <- lam*max(S1_[seq(idx, length(S1_))]) + idx
    }
    v2 <- min(cost)
    i2 <- which(cost==v2)
    if (length(i2)==1){
      r1 <- i2-1
    } else {
      r1 <- max(i2)-1
    }
    lam <- lam + 0.05
    if (itcounter > maxiter){
      break
    }
  }
  
  if (itcounter<=maxiter){
    cost2 <- array(0, c(1, (length(S0)-1)))
    for (idx in seq_len(length(S0)-1)){
      cost2[idx] <- (S0[idx+1]+sqrt(idx * epsilon) * S0[1]/epsilon)/S0[idx]
    }
    v2 <- min(cost2)
    i2 <- which(cost2==v2)
    if (length(i2)==1){
      r2 <- i2
    } else {
      r2 <- max(i2)
    }
    
    if (r1>r2){
      r <- r1
    } else {
      r <- r2
    }
    return(r)
  } else {
    r <- min(nrow(X), ncol(X))
  }
}


# Aux 2 : compute the distortion ------------------------------------------
# @keywords internal
.aux_G <- function(X, m0, r)
{
  z <- rowSums(X^2)/(2*m0*r)
  y <- exp((z-1)^2) - 1
  idxfind <- (z < 1)
  y[idxfind] <- 0
  out <- sum(y)
  return(out)
}

# @keywords internal
.aux_F_t <- function(X, Y, S, M_E, E, m0, rho)
{
  n <- nrow(X)
  r <- ncol(X)
  out1 <- (sum((((X %*% S %*% t(Y)) - M_E)*E)^2))/2
  out2 <- rho*.aux_G(Y,m0,r)
  out3 <- rho*.aux_G(X,m0,r)
  out  <- out1+out2+out3
  return(out)
}


# Aux 3 : compute the gradient --------------------------------------------
# @keywords internal
.aux_Gp <- function(X,m0,r)
{
  z <- rowSums(X^2)/(2*m0*r)
  z <- 2*exp((z-1)^2)/(z-1)
  idxfind <- (z<0)
  z[idxfind] <- 0  
  out <- (X * matrix(z, nrow=nrow(X), ncol=ncol(X), byrow=FALSE))/(m0*r)
}
# @keywords internal
.aux_gradF_t <- function(X, Y, S, M_E, E, m0, rho)
{
  n <- nrow(X)
  r <- ncol(X)
  m <- nrow(Y)
  if (ncol(Y)!=r){
    stop("dimension error from .aux_gradF_t")
  }
  
  XS  <- (X %*% S)
  YS  <- (Y %*% t(S))
  XSY <- (XS %*% t(Y))
  
  Qx <- ((t(X) %*% ((M_E-XSY)*E) %*% YS)/n)
  Qy <- ((t(Y) %*% t((M_E-XSY)*E) %*% XS)/m)
  
  W <- (((XSY-M_E)*E) %*% YS)  + (X %*% Qx) + rho*.aux_Gp(X, m0, r)
  Z <- (t((XSY-M_E)*E) %*% XS) + (Y %*% Qy) + rho*.aux_Gp(Y, m0, r)
  
  resgrad <- list()
  resgrad$W <- W
  resgrad$Z <- Z
  return(resgrad)
  
}


# Aux 4 : Sopt given X and Y ----------------------------------------------
# @keywords internal
.aux_getoptS <- function(X, Y, M_E, E)
{
  n <- nrow(X)
  r <- ncol(X)  
  C <- (t(X) %*% (M_E) %*% Y)
  C <- matrix(as.vector(C))  
  nnrow <- ncol(X)*ncol(Y)
  A <- matrix(NA, nrow=nnrow, ncol=(r^2))
  
  for (i in seq_len(r)){
    for (j in seq_len(r)){
      ind <- ((j-1)*r+i)
      tmp <- t(X) %*% (outer(X[,i], Y[,j])*E) %*% Y      
      A[,ind] <- as.vector(tmp)
    }
  }
  
  S <- solve(A, C)
  out <- matrix(S, nrow=r)
  return(out)
}

# Aux 5 : optimal line search ---------------------------------------------
# @keywords internal
.aux_getoptT <- function(X, W, Y, Z, S, M_E, E, m0, rho)
{
  norm2WZ <- (norm(W, 'f')^2) + (norm(Z, 'f')^2)
  f <- array(0, c(1, 21))
  f[1] <- .aux_F_t(X, Y, S, M_E, E, m0, rho)
  t <- -1e-1
  for (i in seq_len(21)){
    f[i+1] <- .aux_F_t(X+t*W, Y+t*Z, S, M_E, E, m0, rho)
    if ((f[i+1]-f[1]) <= 0.5*t*norm2WZ){
      out <- t
      break
    }
    t <- t/2
  }
  out <- t
  return(t)
}


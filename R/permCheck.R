`permCheck` <- function(object, control = permControl(),
                        make.all = TRUE)
{
    ## if object is numeric or integer and of length 1,
    ## extend the object
    if(length(object) == 1 &&
       (is.integer(object) || is.numeric(object)))
        object <- seq_len(object)
    ## check the number of observations in object
    nobs <- getNumObs(object)
    type <- control$type
    ## if strata, check nobs == length of strata
    ## but beware empty levels
    if(!is.null(control$strata)) {
        tab <- table(control$strata)
        if(!identical(as.integer(nobs), as.integer(sum(tab))))
            stop("Number of observations and length of 'strata' do not match.")
        ## if "grid", check design balanced?
        if((bal <- length(unique(tab))) > 1 && type == "grid")
            stop("Unbalanced 'grid' designs are not supported.")
        ## if grid design, check nrow*ncol is multiple of nobs
        if(type == "grid" &&
           !identical(nobs %% (control$ncol * control$nrow), 0))
            stop("'nrow' * 'ncol' not a multiple of number of observations.")
        ## if constant, check design balanced?
        if(control$constant && bal > 1)
            stop("Unbalanced designs not allowed with 'constant = TRUE'.")
        ## if permuting strata, must be balanced
        if(control$permute.strata && bal > 1)
            stop("Design must be balanced if permuting 'strata'.")
    }
    ##
    if(!is.null(control$all.perms) &&
       !identical(class(control$all.perms), "allPerms"))
        stop("'control$all.perms' must be of class 'allPerms'.")
    ## get number of possible permutations
    num.pos <- numPerms(object, control)
    ## if number of possible perms < minperm turn on complete enumeration
    if(num.pos < control$minperm) {
        control$nperm <- control$maxperm <- num.pos
        control$complete <- TRUE
    }
    ## if complete enumeration, generate all permutations
    if(control$complete && make.all) {
        control$all.perms <- allPerms(nobs, control = control,
                                      max = control$maxperm,
                                      observed = FALSE)
    }
    retval <- list(n = num.pos, control = control)
    class(retval) <- "permCheck"
    retval
}

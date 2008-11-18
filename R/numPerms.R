`numPerms` <- function(object, control = permControl())
{
    ## expand object if a numeric or integer vector of length 1
    if((is.numeric(object) || is.integer(object)) && (length(object) == 1))
         object <- seq_len(object)
    ## number of observations in data
    nobs <- getNumObs(object)
    ## are strata present?
    use.strata <- !is.null(control$strata)
    ## check that when permuting strata or constant within strata,
    ## strata have same number of samples
    if(use.strata) {
        tab.strata <- table(control$strata)
        same.n <- length(unique(tab.strata))
        if((control$permute.strata && same.n > 1) ||
           (control$constant == TRUE && same.n > 1))
            stop("All levels of strata must have same number of samples for chosen scheme")
        if(control$type == "grid" && same.n > 1)
            stop("Unbalanced grid designs are not supported")
    }
    ## generate multiplier for restricted permutations
    if(control$type %in% c("series","grid")) {
        multi <- 2
        if(control$type == "grid" && control$ncol > 2) {
            multi <- 4
        } else {
            if(nobs == 2)
                multi <- 1
        }
    }
    ## calculate number of possible permutations
    num.pos <- if(control$permute.strata) {
        if(control$type == "free")
            exp(lfactorial(length(levels(control$strata))))
        else if(control$type %in% c("series","grid")) {
            if(control$mirror)
                multi * nobs
            else
                nobs
        }
    } else {
        if(control$type == "free") {
            if(use.strata)
                prod(factorial(tab.strata))
            else
                exp(lfactorial(nobs))
        } else if(control$type %in% c("series","grid")) {
            ##multi <- 2
            ##if(control$type == "grid") {
            ##    if(control$ncol == 2)
            ##        multi <- 2
            ##    else
            ##        multi <- 4
            ##} else {
            ##    if(nobs == 2)
            ##        multi <- 1
            ##}
            if(use.strata) {
                if(same.n > 1) {
                    multi <- rep(2, length = length(tab.strata))
                    multi[which(tab.strata == 2)] <- 1
                    if(control$mirror) {
                        prod(multi * tab.strata)
                    } else {
                        prod(tab.strata)
                    }
                } else {
                    if(control$mirror) {
                        if(control$constant)
                            multi * unique(tab.strata)
                        else
                            prod(multi * tab.strata)
                    } else {
                        if(control$constant)
                            unique(tab.strata)
                        else
                            prod(tab.strata)
                    }
                }
            } else {
                if(control$mirror)
                    multi * nobs
                else
                    nobs
            }
        } else {
            stop("Ambiguous permutation type in 'control$type'")
        }
    }
    num.pos
}

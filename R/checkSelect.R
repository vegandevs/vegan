## internal function for checking select arguments in ordination
## plotting functions. Function checks that the 'select' occurs in
## scores (index limits), and takes care that scores are not reordered
## or duplicated. So it returns a set, and if a set is empty it
## returns unmodified scores.
.checkSelect <-
    function(select, scores)
{
    ## select can be a logical vector matching rows of scores, or
    ## vector of numeric indices or vector of rownames of those rows
    ## that are returned. Invalid select entries are silently ignored.
    ##
    ## CHECK: Assumes that scores is either a matrix or a vector and
    ## may fail if scores is a data frame.

    if (is.logical(select)) {
        if (!identical(length(select), NROW(scores))) {
            warning("length of 'select' does not match the number of scores: ignoring 'select'",
                    call. = FALSE)
            return(scores)
        }
        pick <- which(select) # numerical indices
    } else if (is.numeric(select)) {
        N <- if (is.matrix(scores)) NROW(scores) else length(scores)
        pick <- which(seq_len(N) %in% select) # indices within limits
    } else if (is.character(select)) {
        ids <- if (is.matrix(scores)) rownames(scores) else names(scores)
        pick <- which(ids %in% select) # indices in rownames
    } else {
        warning("unknown types in 'select': ignoring 'select'",
                call. = FALSE)
        return(scores)
    }
    if (length(pick) < 1) {
        warning("'select' set is empty: ignoring 'select'", call. = FALSE)
        return(scores)
    }
    if(is.matrix(scores)) {
        scores[pick, , drop = FALSE]
    } else { # vector
        scores[pick]
    }
}

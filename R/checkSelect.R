## internal function for checking select arguments in ordination plotting
## functions
.checkSelect <- function(select, scores) {
    ## check `select` and length of scores match
    if(is.logical(select) &&
                 !isTRUE(all.equal(length(select), NROW(scores)))) {
        warning("length of 'select' does not match the number of scores: ignoring 'select'")
    } else {
        scores <- if(is.matrix(scores)) {
            scores[select, , drop = FALSE]
        } else {
            scores[select]
        }
    }
    scores
}

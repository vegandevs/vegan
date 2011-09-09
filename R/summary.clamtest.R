summary.clamtest <- function(object, ...) {
    structure(c(attr(object, "settings"), 
        list(summary=cbind(Species=table(object$Classes), 
            Proportion=table(object$Classes)/nrow(object)),
        minv=attr(object, "minv"),
        coverage=attr(object, "coverage"))), class="summary.clamtest")
}

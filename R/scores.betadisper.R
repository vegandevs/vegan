`scores.betadisper` <-
    function(x, display = c("sites", "centroids"),
             choices = c(1,2), ...)
{
    display <- match.arg(display, several.ok = TRUE)
    sol <- list()
    if("sites" %in% display)
        sol$sites <- x$vectors[, choices]
    if("centroids" %in% display) {
        if(is.matrix(x$centroids))
            sol$centroids <- x$centroids[, choices, drop = FALSE]
        else
            sol$centroids <- matrix(x$centroids[choices], ncol = length(choices), byrow = TRUE)
    }
    if (length(sol) == 1)
        sol <- sol[[1]]
    sol
}

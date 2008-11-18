`scores.betadisper` <- function(x, display = c("sites", "centroids"),
                                choices = c(1,2), ...)
{
    display <- match.arg(display, several.ok = TRUE)
    #tabula <- c("sites", "centroids")
    #names(tabula) <- c("sites", "centroids")
    #if(length(display) == 1) {
    #    display <- match.arg(display, c("sites", "centroids",
    #                                    "wa", "cn"))
    #    if(display == "sites")
    #        display <- "wa"
    #}
    #take <- tabula[display]
    sol <- list()
    if("sites" %in% display)
        sol$sites <- x$vectors[, choices]
    if("centroids" %in% display)
        sol$centroids <- x$centroids[, choices]
    if (length(sol) == 1) 
        sol <- sol[[1]]
    return(sol)
}

### Written as a helper function to MDSaddpoints, but may be usable
### elsewhere. Can be used to arrange dissmilarities and extract
### partitions of dissimilarity matrix of type

### row: partitions
### 1:
### 2: a
### 3: a a
### 4: a a a
### 5: b b b b
### 6: b b b b c
### 7: b b b b c c

### where the component 'b' is the one we need in MDSaddpoints (and 'a'
### may also be practical).

#' Extract Partitions of Dissimilarities
#'
#' Rectangular partition returns sampling units in \code{pick} against
#' those not in \code{pick} which is identical to distances between
#' observations in two data sets. Symmetric partitions are
#' dissimilarities among observations in \code{}, and are subsets of
#' dissmilarities. The rectangular partitions can be used to specify
#' new observations in \code{\link{MDSaddpoints}}.

#' @examples
#' ## Cross-validation: remove points when performing NMDS and add as
#' ## a new points
#' data(dune)
#' d <- vegdist(dune)
#' ## remove point 3 from ordination
#' mod3 <- metaMDS(dist2xy(d, 3, "xx", reverse = TRUE))
#' ## add point 3 to the result
#' MDSaddpoints(mod3, dist2xy(d, 3))
#'
#'
#' @param dist Input dissimilarities.
#' @param pick Indices (integers) of selected observations. The output
#'     be in the order of the input and will not be reordered by this
#'     argument.
#' @param type \code{"xy"} returns rectangular data of not picked
#'     against picked observations, and \code{"xx"} a subset of
#'     symmetric dissimilarities.
#' @param invert Invert \code{pick}, or drop elements listed.
#'
#' @rdname MDSaddpoints
#' @export
`dist2xy` <-
    function(dist, pick, type = c("xy", "xx"), invert = FALSE)
{
    type <- match.arg(type)
    if (!inherits(dist, "dist"))
        stop("'dist' must be a dissimilarity object of class 'dist'")
    att <- attributes(dist)
    dist <- as.matrix(dist)
    n <- nrow(dist)
    ## make up the selection vector
    k <- logical(n)
    k[pick] <- TRUE
    if (invert)
        k <- !k
    ## make compartments
    dist <- if (type == "xy")
                dist[k, !k, drop = FALSE]
            else
                as.dist(dist[k, k, drop = FALSE])
    attr(dist, "call") <- match.call()
    attr(dist, "method") <- att$method
    attr(dist, "maxdist") <- att$maxdist
    dist
}


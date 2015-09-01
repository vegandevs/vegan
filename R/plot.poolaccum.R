`plot.poolaccum` <-
    function(x, alpha = 0.05, type = c("l","g"), ...)
{
    m <- summary(x, alpha = alpha, ...)
    n <- nrow(m[[1]])
    Size <- as.vector(sapply(m, function(x) c(x[,1], x[,1], rev(x[,1]))))
    Richness <- as.vector(sapply(m, function(x) c(x[,2], x[,3], rev(x[,4]))) )
    indnames <- as.character(sapply(m, function(x) colnames(x[,2, drop=FALSE])))
    Index <- factor(rep(indnames, each = 3 * n), levels = indnames)
    lintype <- rep(c(rep("aver", n), rep("envel", 2*n)), length=length(Size))
    xyplot(Richness ~  Size | Index, as.table = TRUE, groups = lintype,
           type = type, ...)
}

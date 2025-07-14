`distconnected` <-
    function(dis, toolong = 1, trace = TRUE)
{
    n <- attr(dis, "Size")
    out <- .C(stepabyss, dis = as.double(dis), n = as.integer(n),
              toolong = as.double(toolong), val = integer(n),
              NAOK = TRUE, PACKAGE = "vegan")$val
    if (trace) {
        cat("Connectivity of distance matrix with threshold dissimilarity",
            toolong,"\n")
        n <- length(unique(out))
        if (n == 1)
            cat("Data are connected\n")
        else {
            cat("Data are disconnected:", n, "groups\n")
            print(table(out, dnn="Groups sizes"))
        }
    }
    out
}

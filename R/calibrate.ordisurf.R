`calibrate.ordisurf` <-
    function(object, newdata, ...)
{
    if (missing(newdata))
        fit <- predict(object, type = "response", ...)
    else {
        ## Got only a vector of two coordinates
        if (is.vector(newdata) && length(newdata) == 2) 
            newdata = data.frame(x1 = newdata[1], x2 = newdata[2])
        ## Got a matrix or a data frme
        else{
            if (NCOL(newdata) < 2) 
                stop("needs a matrix or a data frame with two columns")
            newdata <- data.frame(x1 = newdata[,1], x2 = newdata[,2])
        }
        fit <- predict(object, newdata = newdata, type = "response", ...)
    }
    fit
}

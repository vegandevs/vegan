### predict method for radline, radfit & radfit.frame

### All functions take 'newdata' argument which need not be integer:
### the functions can interpolate, but not necessarily extrapolate, or
### the extrapolations may be NaN.

`predict.radline`  <-
    function(object, newdata, ...)
{
    if (missing(newdata))
        x <- seq_along(object$y)
    else
        x <- drop(as.matrix(newdata))
    nobs <- length(object$y)
    p <- coef(object)
    switch(object$model,
           ## linear interpolation, no extrapolation
           `Brokenstick` = approx(seq_len(nobs), object$fitted.values, x, ...)$y,
           `Preemption` = exp(log(sum(object$y)) + log(p) + log(1 - p)*(x-1)),
           ## NaN when rank outside proportional rank 0...1 
           `Log-Normal` = {
               slope <- diff(range(ppoints(nobs)))/(nobs-1)
               intcpt <- 0.5 - slope * (nobs + 1) / 2
               xnorm <- -qnorm(intcpt + slope * x)
               exp(p[1] + p[2]*xnorm)
           },
           `Zipf` = exp(log(sum(object$y)) + log(p[1]) + p[2]*log(x)),
           `Zipf-Mandelbrot` = exp(log(sum(object$y)) + log(p[1]) +
           p[2]*log(x + p[3]))
           )
}

`predict.radfit`<-
    function(object, newdata, ...)
{
    if (missing(newdata))
        sapply(names(object$models), function(x, ...)
               predict(object$models[[x]], ...))
    else
        sapply(names(object$models), function(x, ...)
               predict(object$models[[x]], newdata, ...))
}

`predict.radfit.frame`  <-
    function(object, newdata, ...)
{
    if(missing(newdata))
        lapply(object, predict, ...)
    else
        lapply(object, predict, newdata = newdata)
}

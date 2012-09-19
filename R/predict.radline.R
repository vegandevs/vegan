### predict method for radline, radfit & radfit.frame

### All functions take 'newdata' argument which need not be integer:
### the functions can interpolate, but not necessarily extrapolate, or
### the extrapolations may be NaN.

`predict.radline`  <-
    function(object, newdata, total, ...)
{
    ## newdata can be ranks
    if (missing(newdata))
        x <- seq_along(object$y)
    else
        x <- drop(as.matrix(newdata))
    ## total number of individuals in the community
    if (missing(total))
        total <- sum(object$y)
    ## adjustment for chagned total in call
    adj <- total/sum(object$y)
    nobs <- length(object$y)
    p <- coef(object)
    switch(object$model,
           ## linear interpolation, no extrapolation
           `Brokenstick` = approx(seq_len(nobs),
           object$fitted.values, x, ...)$y * adj,
           `Preemption` = exp(log(total) + log(p) + log(1 - p)*(x-1)),
           ## NaN when rank outside proportional rank 0...1 
           `Log-Normal` = {
               slope <- diff(range(ppoints(nobs)))/(nobs-1)
               intcpt <- 0.5 - slope * (nobs + 1) / 2
               xnorm <- -qnorm(intcpt + slope * x)
               exp(p[1] + p[2]*xnorm)*adj
           },
           `Zipf` = exp(log(total) + log(p[1]) + p[2]*log(x)),
           `Zipf-Mandelbrot` = exp(log(total) + log(p[1]) +
           p[2]*log(x + p[3]))
           )
}

`predict.radfit`<-
    function(object, ...)
{
    sapply(object$models, predict, ...)
}

`predict.radfit.frame`  <-
    function(object, ...)
{
    lapply(object, predict, ...)
}

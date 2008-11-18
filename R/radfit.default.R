"radfit.default" <-
    function (x, ...) 
{
    x <- as.rad(x)
    NU <- rad.null(x, ...)
    PE <- rad.preempt(x, ...)
    ##BS <- rad.brokenstick(x, ...)
    LN <- rad.lognormal(x, ...)
    ZP <- rad.zipf(x, ...)
    ZM <- rad.zipfbrot(x, ...)
    out <- list(y = x, family = PE$family)
    models <- list(Null = NU, Preemption = PE, Lognormal = LN,
                   Zipf = ZP, Mandelbrot = ZM)
    out$models <- models
    class(out) <- "radfit"
    out
}

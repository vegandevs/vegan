`panel.ordi3d` <-
    function(x, y, z, aspect, ...)
{
    dx <- diff(range(x))
    dy <- diff(range(y))
    dz <- diff(range(z))
    aspect <- c(dy/dx, dz/dx)
    panel.cloud(x = x, y = y, z = z, aspect = aspect, ...)
}

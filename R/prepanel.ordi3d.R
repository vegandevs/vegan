`prepanel.ordi3d` <-
    function(xlim = xlim, ylim = ylim, zlim = zlim, aspect = c(1,1),  ...)
{
    aspect = c(diff(ylim)/diff(xlim), diff(zlim)/diff(xlim))
    prepanel.default.cloud(xlim = xlim, ylim = ylim, zlim = zlim, aspect = aspect, ...)
}

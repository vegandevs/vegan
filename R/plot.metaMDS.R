"plot.metaMDS" <-
    function (x,  display = c("sites","species"), choices = c(1, 2), type = "p",
              shrink = FALSE, ...) 
{
    if (length(display) == 1)
        display <- match.arg(display)
    if (shrink) {
        x$species <- scores(x, display = "species", shrink = shrink,
                            choices=1:x$dims)
    }
    ordiplot(x, choices = choices, type = type, display = display, ...)
}

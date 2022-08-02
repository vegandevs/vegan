`plot.metaMDS`  <-
    function (x,  display = c("sites","species"), choices = c(1, 2), type = "p",
              shrink = FALSE, ...)
{
    x <- scores(x, display = display, choices = choices, shrink = shrink)
    ordiplot(x, display = display, type = type, ...)
}

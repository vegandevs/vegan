"eigengrad" <-
function (x, w) 
{
    attr(wascores(x, w, expand=TRUE), "shrinkage")
}

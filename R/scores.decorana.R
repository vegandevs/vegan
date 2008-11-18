"scores.decorana" <-
function (x, display=c("sites","species"), choices = 1:4, origin=TRUE, ...) 
{
   display <- match.arg(display)
   if(display == "sites")
      X <- x$rproj
   else if(display == "species")
      X <- x$cproj
   if (origin)
      X <- sweep(X, 2, x$origin, "-")
   X <- X[,choices]
   X 
}

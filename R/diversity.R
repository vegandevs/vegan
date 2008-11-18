"diversity" <-
  function (x, index = "shannon", MARGIN = 1, base = exp(1)) 
{
  x <- as.matrix(x)
  INDICES <- c("shannon", "simpson", "invsimpson")
  index <- match.arg(index, INDICES)
  total <- apply(x, MARGIN, sum)
  x <- sweep(x, MARGIN, total, "/")
  if (index == "shannon")
    x <- -x * log(x, base)
  else
    x <- x^2
  H <- apply(x, MARGIN, sum, na.rm = TRUE)
  if (index == "simpson") 
    H <- 1 - H
  else if (index == "invsimpson") 
    H <- 1/H
  return(H)
}



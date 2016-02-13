"initMDS" <-
  function(x, k=2)
{
  nr <- attr(x, "Size")
  res <- runif(nr*k)
  dim(res) <- c(nr,k)
  res
}

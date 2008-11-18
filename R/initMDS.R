"initMDS" <-
  function(x, k=2)
{
  nr <- attributes(x)$Size
  res <- runif(nr*k)
  dim(res) <- c(nr,k)
  res
}

"wisconsin" <-
  function (x) 
{
  x <- decostand(x, "max", 2)
  x <- decostand(x, "tot", 1)
  x 
}

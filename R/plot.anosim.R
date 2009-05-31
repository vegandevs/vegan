"plot.anosim" <-
function (x, title=NULL, ...) 
{
   boxplot(x$dis.rank ~ x$class.vec, notch=TRUE, varwidth=TRUE, 
           ...)
   title(title)
   if (x$permutations) {
     pval <- format.pval(x$signif)
   } else {
     pval <- "not assessed"
   }
   mtext(paste("R = ", round(x$statistic, 3), ", ",
               "P = ", pval ), 3)
   invisible()
}

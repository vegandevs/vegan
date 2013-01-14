`summary.anosim` <-
function (object, ...) 
{
   print(object)
   if (object$permutations) {
     out <- quantile(object$perm, c(0.9, 0.95, 0.975, 0.99))
     cat("Upper quantiles of permutations (null model):\n")
     print(out, digits=3)
   }
   cat("\n")
   tmp <- tapply(object$dis.rank, object$class.vec, quantile)
   out <- matrix(NA, length(tmp), 5)
   for (i in 1:length(tmp)) out[i,] <- tmp[[i]]
   rownames(out) <- names(tmp)
   colnames(out) <- names(tmp$Between)
   out <- cbind(out, N = table(object$class.vec))
   cat("Dissimilarity ranks between and within classes:\n")
   print(out)
   cat("\n") 
   invisible()
}

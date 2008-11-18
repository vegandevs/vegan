`summary.taxondive` <-
function (object, ...) 
{ 
    z <- (object$Dplus - object$EDplus)/object$sd.Dplus
    pval <- 2*pnorm(-abs(z))
    out <- cbind(object$D, object$Dstar, object$Dplus, object$sd.Dplus, 
                 z, pval)
    out <- rbind(out, "Expected"=c(object$ED, object$EDstar, object$EDplus, NA, NA, NA))
    colnames(out) <- c("Delta", "Delta*", "Delta+", "sd(Delta+)", 
                       "z(Delta+)", "Pr(>|z|)")
    class(out) <- "summary.taxondive"
    out
}


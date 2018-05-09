`anova.betadisper` <- function(object, ...)
{
    model.dat <- with(object, data.frame(Distances = distances,
                                         Groups = group))
    n.grps <- with(model.dat, length(unique(as.numeric(Groups))))
    if(n.grps < 2)
        stop("anova() only applicable to two or more groups")
    anova(lm(Distances ~ Groups, data = model.dat))
}

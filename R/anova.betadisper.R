`anova.betadisper` <- function(object, ...)
{
    model.dat <- with(object, data.frame(Distances = distances,
                                         Groups = group))
    anova(lm(Distances ~ Groups, data = model.dat))
}

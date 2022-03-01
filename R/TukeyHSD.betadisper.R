`TukeyHSD.betadisper` <- function(x, which = "group", ordered = FALSE,
                                  conf.level = 0.95, ...) {
    df <- data.frame(distances = x$distances, group = x$group)
    mod.aov <- aov(distances ~ group, data = df)
    TukeyHSD(mod.aov, which = which, ordered = ordered,
             conf.level = conf.level, ...)
}

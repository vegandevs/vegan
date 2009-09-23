\name{ordilabel}
\alias{ordilabel}

\title{Add Text on Non-transparent Label to an Ordination Plot }
\description{
 Function \code{ordilabel} is similar to \code{\link{text}}, but the text is on an
 opaque label. This can help in crowded ordination plots: you still cannot see
 all text labels, but at least the uppermost are readable. Argument \code{priority}
 helps to make the most important labels visible.
}
\usage{
ordilabel(x, display, labels, choices = c(1, 2), priority, cex = 0.8,
    fill = "white", border = NULL, ...)
}

\arguments{
  \item{x}{An ordination object an any object known to \code{\link{scores}}. }
  \item{display}{Kind of scores displayed (passed to \code{\link{scores}}). }
  \item{labels}{Optional text used in plots. If this is not given, the text is found from the
       ordination object.}
  \item{choices}{Axes shown (passed to \code{\link{scores}}). }
  \item{priority}{Vector of the same length as the number of labels. The items with
       high priority will be plotted uppermost. }
  \item{cex}{ Character expansion for the text (passed to \code{\link{text}}). }
  \item{fill}{ Background colour of the labels (the \code{col} argument of
        \code{\link{polygon}}).}
  \item{border}{The colour and visibility of the border of the label as defined in
         \code{\link{polygon}}).}
  \item{\dots}{Other arguments (passed to \code{\link{text}}). }
}
\details{
 The function may be useful with crowded ordination plots, in particular together with
  argument \code{priority}. You will not see all text labels, but at least some are
  readable. Other alternatives to crowded plots are
  \code{\link{identify.ordiplot}}, \code{\link{orditorp}} and \code{\link{orditkplot}}. 
}

\author{ Jari Oksanen }

\seealso{ \code{\link{scores}}, \code{\link{polygon}}, \code{\link{text}}. The function is
 modelled after \code{\link[ade4]{s.label}} in \pkg{ade4} package.}
\examples{
data(dune)
ord <- cca(dune)
plot(ord, type = "n")
ordilabel(ord, dis="sites", cex=1.2, font=3, fill="hotpink", col="blue")
## You may prefer separate plots, but here species as well
ordilabel(ord, dis="sp", font=2, priority=colSums(dune))
}
\keyword{ aplot }


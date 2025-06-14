\name{ordilabel}
\alias{ordilabel}

\title{Add Text on Non-transparent Label to an Ordination Plot. }

\description{ Function \code{ordilabel} is similar to
  \code{\link{text}}, but the text is on an opaque label. This can help
  in crowded ordination plots: you still cannot see all text labels, but
  at least the uppermost ones are readable. Argument \code{priority} helps to
  make the most important labels visible.  Function can be used in pipe
  after ordination \code{plot} or \code{\link{ordiplot}} command. }

\usage{
ordilabel(x, display, labels, choices = c(1, 2), priority, select,
    cex = 0.8, fill = "white", border = NULL, col = NULL, xpd = TRUE, ...)
}

\arguments{
  \item{x}{An ordination object an any object known to
    \code{\link{scores}}. }
  \item{display}{Kind of scores displayed (passed to
    \code{\link{scores}}). }
  
  \item{labels}{Optional text used in plots instead of the default. If
    \code{select} is given, the labels are given only to selected items
    in the order they occur in the scores.}

  \item{choices}{Axes shown (passed to \code{\link{scores}}). }

  \item{priority}{Vector of the same length as the number of scores or
    selected items. The items with high priority will be plotted
    uppermost.}

  \item{select}{Items to be displayed. This can either be a logical
    vector which is \code{TRUE} for displayed items or a vector of
    indices of displayed items.}

  \item{cex}{ Character expansion for the text (passed to
    \code{\link{text}}). }
  \item{fill}{ Background colour of the labels (the \code{col} argument
    of \code{\link{polygon}}).}
  \item{border}{The colour and visibility of the border of the label as
    defined in \code{\link{polygon}}. The default is to use text colour
    \code{col}.}
  \item{col}{Text colour.}
  \item{xpd}{Draw labels also outside the plot region.}
  \item{\dots}{Other arguments (passed to \code{\link{text}}). }
}

\details{ The function may be useful with crowded ordination plots, in
  particular together with argument \code{priority}. You will not see
  all text labels, but at least some are readable. Function can be used
  as a part of a pipe (\code{|>}) in place of \code{text} after an
  ordination \code{plot} command (see Examples).

  Other alternatives for cluttered plots are
  \code{\link{identify.ordiplot}}, \code{\link{orditorp}},
  \code{\link{ordipointlabel}}, and \code{orditkplot} (\CRANpkg{vegan3d}
  package).
}

\author{ Jari Oksanen }

\seealso{\code{\link{plot.cca}} and \code{\link{text.ordiplot}} that
  can use the function with argument \code{bg}.}

\examples{
data(dune)
ord <- cca(dune)
plot(ord, type = "n")
## add text
ordilabel(ord, dis="sites", cex=1.2, font=3, fill="hotpink", col="blue")
## You may prefer separate plots, but here species as well
ordilabel(ord, dis="sp", font=2, priority=colSums(dune))
## use in a pipe
plot(ord, type = "n") |>
   ordilabel("spec", font = 3, priority = colSums(dune)) |>
   points("sites", pch=21, bg = "yellow", col = "blue")
}
\keyword{ aplot }


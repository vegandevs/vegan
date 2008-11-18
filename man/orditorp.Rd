\name{orditorp}
\alias{orditorp}

\title{ Add Text or Points to Ordination Plots }
\description{
  The function adds \code{\link{text}} or \code{\link{points}} to
  ordination plots.  Text will be used if this can be done without
  overwriting other text labels, and points will be used otherwise. The
  function can help in reducing clutter in ordination graphics, but
  manual editing may still be necessary.
}
\usage{
orditorp(x, display, labels, choices = c(1, 2), priority,
    cex = 0.7, pcex, col = par("col"), pcol,
    pch = par("pch"), air = 1, ...)
}

\arguments{
  \item{x}{A result object from ordination or an \code{\link{ordiplot}}
    result. }
  \item{display}{Items to be displayed in the plot.  Only one
    alternative is allowed. Typically this is \code{"sites"} or
    \code{"species"}.  }
  \item{labels}{ Optional text used for labels. Row names will be used if
  this is missing. }
  \item{choices}{ Axes shown.}
  \item{priority}{ Text will be used for items with higher priority
    if labels overlap.  This should be vector of the same
    length as the number of items plotted.}
  \item{cex, pcex}{Text and point sizes, see \code{\link{plot.default}}..}
  \item{col, pcol}{Text and point colours, see \code{\link{plot.default}}.}
  \item{pch}{Plotting character, see \code{\link{points}}.}
  \item{air}{Amount of empty space between text labels. Values <1 allow
    overlapping text.}
  \item{\dots}{ Other arguments to \code{\link{scores}} (and its various
    methods), \code{\link{text}} and \code{\link{points}}.}
}
\details{
  Function \code{orditorp} will add either text or points to an existing
  plot. The items with high \code{priority} will be added first
  and \code{\link{text}} will be used if this can be done without
  overwriting previous labels,and \code{\link{points}} will be used
  otherwise. If \code{priority} is missing, labels will be added from the
  outskirts to the centre. Function \code{orditorp} can be used
  with most ordination results, or plotting results from
  \code{\link{ordiplot}} or ordination plot functions
  (\code{\link{plot.cca}}, \code{\link{plot.decorana}},
  \code{\link{plot.metaMDS}}).

  Arguments can be passed to the relevant \code{\link{scores}} method
  for the ordination object (\code{x}) being drawn. See the relevant
  \code{\link{scores}} help page for arguments that can be used.
}
\value{
  The function returns invisibly a logical vector where \code{TRUE}
  means that item was labelled with text and \code{FALSE} means that it
  was marked with a point. The returned vector can be used as the
  \code{select} argument in ordination \code{text} and \code{points}
  functions.
}

\author{ Jari Oksanen }


\examples{
## A cluttered ordination plot :
data(BCI)
mod <- cca(BCI)
plot(mod, dis="sp", type="t")
# Now with orditorp and abbreviated species names
cnam <- make.cepnames(names(BCI))
plot(mod, dis="sp", type="n")
stems <- colSums(BCI)
orditorp(mod, "sp", label = cnam, priority=stems, pch="+", pcol="grey")
}
\keyword{ aplot }
\keyword{ hplot }

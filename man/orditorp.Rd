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
    select, cex = 0.7, pcex, col = par("col"), pcol,
    pch = par("pch"), air = 1, ...)
}

\arguments{
  \item{x}{A result object from ordination or an \code{\link{ordiplot}}
    result. If the function is used in \code{\link{ordiplot}} pipe, this
    should be missing and first argument be \code{display}.}
  \item{display}{Items to be displayed in the plot.  Only one
    alternative is allowed. Typically this is \code{"sites"} or
    \code{"species"}.  }
  \item{labels}{ Optional text used for labels. Row names of scores will
    be used if this is missing. If \code{select} is used, labels are
    given only selected items in the order they occur in the scores. }
  \item{choices}{ Axes shown.}
  \item{priority}{ Text will be used for items with higher priority
    if labels overlap.  This should be vector of the same
    length as the number of items plotted or number of scores.}
  \item{select}{Items to be displayed.  This can either be a logical
    vector which is \code{TRUE} for displayed items or a vector of
    indices, or \code{labels} and if these are missing, row names of
    scores. If a logical vector is used, it must have the same length as
    the scores plotted.} 
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
  \code{\link{plot.metaMDS}}). Function can also be used in a pipe
  (\code{|>}) where the first command is a \pkg{vegan} ordination
  \code{plot} command or \code{\link{ordiplot}}.

  Arguments can be passed to the relevant \code{\link{scores}} method
  for the ordination object (\code{x}) being drawn. See the relevant
  \code{\link{scores}} help page for arguments that can be used.
}
\value{
  The function returns invisibly the The function returns invisibly a
  logical vector where \code{TRUE} means that item was labelled with
  text and \code{FALSE} means that it was marked with a point. If
  function is used in an \code{\link{ordiplot}} pipe, it will return the
  input \code{ordiplot} object, but amend the plotted scores with this
  vector as attribute \code{"orditorp"}. The returned vector can be used
  as the \code{select} argument in ordination \code{text} and
  \code{points} functions.
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
orditorp(mod, "sp", labels = cnam, priority=stems, pch="+", pcol="grey")

## show select in action
set.seed(1)
take <- sample(ncol(BCI), 50)
plot(mod, dis="sp", type="n")
stems <- colSums(BCI)
## only selected items are labelled, and the labels must be in the some
## order as in the scores
orditorp(mod, "sp", labels = cnam[sort(take)], priority=stems, select = take,
         pch="+", pcol="grey")
\dontshow{
## example(orditorp) should not set random seed in the user session
rm(.Random.seed)
}
}
\keyword{ aplot }
\keyword{ hplot }

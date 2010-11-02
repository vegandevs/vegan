\name{plot.cca}
\alias{plot.cca}
\alias{text.cca}
\alias{points.cca}
\alias{scores.cca}
\alias{scores.rda}
\alias{summary.cca}
\alias{print.summary.cca}
\alias{ade2vegancca}
\alias{head.summary.cca}
\alias{tail.summary.cca}

\title{Plot or Extract Results of Constrained Correspondence Analysis
  or Redundancy Analysis}
\description{
  Functions to plot or extract results of constrained correspondence analysis
  (\code{\link{cca}}), redundancy analysis (\code{\link{rda}}) or
  constrained analysis of principal coordinates (\code{\link{capscale}}).
}
\usage{
\method{plot}{cca}(x, choices = c(1, 2), display = c("sp", "wa", "cn"),
         scaling = 2, type, xlim, ylim, const, ...)
\method{text}{cca}(x, display = "sites", labels, choices = c(1, 2), scaling = 2,
    arrow.mul, head.arrow = 0.05, select, const, ...)
\method{points}{cca}(x, display = "sites", choices = c(1, 2), scaling = 2,
    arrow.mul, head.arrow = 0.05, select, const, ...)
\method{scores}{cca}(x, choices=c(1,2), display=c("sp","wa","cn"), scaling=2, ...)
\method{scores}{rda}(x, choices=c(1,2), display=c("sp","wa","cn"), scaling=2, 
    const, ...)
\method{summary}{cca}(object, scaling = 2, axes = 6, display = c("sp", "wa", 
    "lc", "bp", "cn"), digits = max(3, getOption("digits") - 3), ...)
\method{print}{summary.cca}(x, digits = x$digits, head = NA, tail = head, ...)
\method{head}{summary.cca}(x, n = 6, tail = 0, ...)
\method{tail}{summary.cca}(x, n = 6, head = 0, ...)
}

\arguments{
  \item{x, object}{A \code{cca} result object.}
    \item{choices}{Axes shown.}
  \item{display}{Scores shown.  These must some of the alternatives
    \code{species} or \code{sp} for species scores, \code{sites} or
    \code{wa} for site scores, \code{lc} for linear constraints or ``LC
    scores'', or \code{bp} for biplot arrows or \code{cn} for centroids
    of factor constraints instead of an arrow.}
   \item{scaling}{Scaling for species and site scores. Either species
    (\code{2}) or site (\code{1}) scores are scaled by eigenvalues, and
    the other set of scores is left unscaled, or with \code{3} both are
    scaled symmetrically by square root of eigenvalues.  Corresponding
    negative values can be used in \code{cca} to additionally multiply
    results with \eqn{\sqrt(1/(1-\lambda))}.  This scaling is know as
    Hill scaling (although it has nothing to do
    with Hill's rescaling of \code{\link{decorana}}). With corresponding
    negative values in\code{rda},
    species scores are divided by standard deviation of each species and
    multiplied with an equalizing constant.
    Unscaled raw scores stored in the result can be accessed with
    \code{scaling = 0}. 
  }
  \item{type}{Type of plot: partial match to \code{text}
    for text labels, \code{points} for points, and \code{none} for
    setting frames only.  If omitted, \code{text} is selected for
    smaller data sets, and \code{points} for larger.}
  \item{xlim, ylim}{the x and y limits (min,max) of the plot.}
  \item{labels}{Optional text to be used instead of row names.}
  \item{arrow.mul}{Factor to expand arrows in the graph.  Arrows will be
    scaled automatically to fit the graph if this is missing.}
  \item{head.arrow}{Default length of arrow heads.}
  \item{select}{Items to be displayed.  This can either be a logical
    vector which is \code{TRUE} for displayed items or a vector of indices
    of displayed items.}
  \item{const}{General scaling constant to \code{rda} scores. The
    default is to use a constant that gives biplot scores, that is,
    scores that approximate original data (see \code{\link{vignette}}
    \samp{decision-vegan.pdf} with \code{\link{vegandocs}} for details
    and discussion). If \code{const} is a vector of two items, the first
    is used for species, and the second item for site scores.}
  \item{axes}{Number of axes in summaries.}
  \item{digits}{Number of digits in output.}
  \item{n, head, tail}{Number of rows printed from the head and tail of
    species and site scores.  Default \code{NA} prints all.}
  \item{...}{Parameters passed to other functions.}
}

\details{
  Same \code{plot} function will be used for \code{\link{cca}} and
  \code{\link{rda}}. This produces a quick, standard plot with current
  \code{scaling}.

  The \code{plot} function sets colours (\code{col}), plotting
  characters (\code{pch}) and character sizes (\code{cex}) to
  certain standard values. For a fuller control of produced plot, it is
  best to call \code{plot} with \code{type="none"} first, and then add
  each plotting item separately using \code{text.cca} or
  \code{points.cca} functions. These use the default settings of standard
  \code{\link{text}} and \code{\link{points}} functions and accept all
  their parameters, allowing  a full user control of produced plots.

  Environmental variables receive a special treatment. With
  \code{display="bp"}, arrows will be drawn. These are labelled with
  \code{text} and unlabelled with \code{points}. The basic \code{plot}
  function uses a simple (but not very clever) heuristics for adjusting
  arrow lengths to plots, but the user can give the expansion factor in
  \code{mul.arrow}. With \code{display="cn"} the centroids of levels of
  \code{\link{factor}} variables are displayed (these are available only if there were
  factors and a formula interface was used in \code{\link{cca}} or
  \code{\link{rda}}). With this option continuous
  variables still are presented as arrows and ordered factors as arrows
  and centroids. 

  If you want to have still a better control of plots, it is better to
  produce them using primitive \code{plot} commands. Function
  \code{scores} helps in extracting the 
  needed components with the selected \code{scaling}.

  Function \code{summary} lists all scores and the output can be very
  long.  You can suppress scores by setting \code{axes = 0} or
  \code{display = NA} or \code{display = NULL}. You can display some
  first or last (or both) rows of scores by using \code{head} or
  \code{tail} or explicit \code{print} command for the \code{summary}.
  
 Palmer (1993) suggested using linear constraints
  (``LC scores'') in ordination diagrams, because these gave better
  results in simulations and site scores (``WA scores'') are a step from
  constrained to unconstrained analysis.  However, McCune (1997) showed
  that noisy environmental variables (and all environmental
  measurements are noisy) destroy ``LC scores'' whereas ``WA scores''
  were little affected.  Therefore the \code{plot} function uses site
  scores (``WA scores'') as the default. This is consistent with the
  usage in statistics and other functions in \R
  (\code{\link[MASS]{lda}}, \code{\link{cancor}}).
}
\value{
  The \code{plot} function returns invisibly a plotting structure which
  can be used by function \code{\link{identify.ordiplot}} to identify
  the points or other functions in the \code{\link{ordiplot}} family. 
}

\note{Package \pkg{ade4} has function \code{\link[ade4]{cca}} which
  returns constrained correspondence analysis of the same class as the
  \pkg{vegan} function. If you have results of \pkg{ade4} in your
  working environment, \pkg{vegan} functions may try to handle them and
  fail with cryptic error messages. However, there is a simple utility
  function \code{ade2vegancca} which tries to translate \pkg{ade4}
  \code{cca} results to \pkg{vegan} \code{cca} results so that some
  \pkg{vegan} functions may work partially with \pkg{ade4} objects
  (with a warning).} 

\author{Jari Oksanen }

\seealso{\code{\link{cca}}, \code{\link{rda}} and \code{\link{capscale}}
  for getting something
  to plot, \code{\link{ordiplot}} for an alternative plotting routine
  and more support functions, and \code{\link{text}},
  \code{\link{points}} and \code{\link{arrows}} for the basic routines.  }

\examples{
data(dune)
data(dune.env)
mod <- cca(dune ~ A1 + Moisture + Management, dune.env)
plot(mod, type="n")
text(mod, dis="cn")
points(mod, pch=21, col="red", bg="yellow", cex=1.2)
text(mod, "species", col="blue", cex=0.8)
## Limited output of 'summary'
head(summary(mod), tail=2)
## Read description of scaling in RDA in vegan:
\dontrun{vegandocs("decision")}}
\keyword{hplot}
\keyword{aplot}

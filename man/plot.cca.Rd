\encoding{UTF-8}
\name{plot.cca}
\alias{plot.cca}
\alias{text.cca}
\alias{points.cca}
\alias{scores.cca}
\alias{scores.rda}
\alias{summary.cca}
\alias{print.summary.cca}
\alias{labels.cca}

\title{Plot or Extract Results of Constrained Correspondence Analysis
  or Redundancy Analysis}
\description{
  Functions to plot or extract results of constrained correspondence analysis
  (\code{\link{cca}}), redundancy analysis (\code{\link{rda}}), distance-based
  redundancy analysis (\code{\link{dbrda}}) or
  constrained analysis of principal coordinates (\code{\link{capscale}}).
}
\usage{
\method{plot}{cca}(x, choices = c(1, 2), display = c("sp", "wa", "cn"),
     scaling = "species", type, xlim, ylim, const,
     correlation = FALSE, hill = FALSE, spe.par = list(), sit.par =
     list(), con.par = list(), bip.par = list(), cen.par = list(),
     reg.par = list(), ...)
\method{text}{cca}(x, display = "sites", labels, choices = c(1, 2),
     scaling = "species", arrow.mul, head.arrow = 0.05, select, const,
     axis.bp = FALSE, correlation = FALSE, hill = FALSE, ...)
\method{points}{cca}(x, display = "sites", choices = c(1, 2),
       scaling = "species", arrow.mul, head.arrow = 0.05, select, const,
       axis.bp = FALSE, correlation = FALSE, hill = FALSE, ...)
\method{scores}{cca}(x, choices = c(1,2), display = "all",
       scaling = "species", hill = FALSE, tidy = FALSE, droplist = TRUE,
       ...)
\method{scores}{rda}(x, choices = c(1,2), display = "all",
       scaling = "species", const, correlation = FALSE, tidy = FALSE,
       droplist = TRUE, ...)
\method{summary}{cca}(object, digits = max(3, getOption("digits") - 3), ...)
\method{labels}{cca}(object, display, ...)
}

\arguments{
  \item{x, object}{A \code{cca} result object.}
    \item{choices}{Axes shown.}
  \item{display}{Scores shown.  These must include some of the
    alternatives \code{"species"} or \code{"sp"} for species scores,
    \code{sites} or \code{"wa"} for site scores, \code{"lc"} for linear
    constraints or LC scores, or \code{"bp"} for biplot arrows or
    \code{"cn"} for centroids of factor constraints instead of an arrow,
    and \code{"reg"} for regression coefficients (a.k.a. canonical
    coefficients). The alternative \code{"all"} selects all available
    scores.}
  \item{scaling}{Scaling for species and site scores. Either species
    (\code{2}) or site (\code{1}) scores are scaled by eigenvalues, and
    the other set of scores is left unscaled, or with \code{3} both are
    scaled symmetrically by square root of eigenvalues.  Corresponding
    negative values can be used in \code{cca} to additionally multiply
    results with \eqn{\sqrt(1/(1-\lambda))}.  This scaling is know as Hill
    scaling (although it has nothing to do with Hill's rescaling of
    \code{\link{decorana}}). With corresponding negative values
    in \code{rda}, species scores are divided by standard deviation of each
    species and multiplied with an equalizing constant. Unscaled raw
    scores stored in the result can be accessed with \code{scaling = 0}.

    The type of scores can also be specified as one of \code{"none"},
    \code{"sites"}, \code{"species"}, or \code{"symmetric"}, which
    correspond to the values \code{0}, \code{1}, \code{2}, and \code{3}
    respectively. Arguments \code{correlation} and \code{hill} in
    \code{scores.rda} and \code{scores.cca} respectively can be used in
    combination with these character descriptions to get the
    corresponding negative value.
  }
  \item{correlation, hill}{logical; if \code{scaling} is a character
    description of the scaling type, \code{correlation} or \code{hill}
    are used to select the corresponding negative scaling type; either
    correlation-like scores or Hill's scaling for PCA/RDA and CA/CCA
    respectively. See argument \code{scaling} for details.}
  
  \item{spe.par, sit.par, con.par, bip.par, cen.par, reg.par}{
    Lists of graphical parameters for species, sites, constraints (lc
    scores), biplot and text, centroids and regression. These take
    precedence over globally set parameters and defaults.} 
  
  \item{tidy}{Return scores that are compatible with
    \CRANpkg{ggplot2}: all scores are in a single \code{data.frame},
    score type is identified by factor variable \code{score}, the
    names by variable \code{label}, and weights (in CCA) are in
    variable \code{weight}. The possible values of \code{score} are
    \code{species}, \code{sites} (for WA scores), \code{constraints}
    (LC scores for sites calculated directly from the constraining
    variables), \code{biplot} (for biplot arrows), \code{centroids}
    (for levels of factor variables), \code{factorbiplot} (biplot
    arrows that model centroids), \code{regression} (for regression
    coefficients to find LC scores from constraints). These scores
    cannot be used with conventional \code{plot}, but they are
    directly suitable to be used with the \pkg{ggplot2} package.}

  \item{type}{Type of plot: partial match to \code{text}
    for text labels, \code{points} for points, and \code{none} for
    setting frames only.  If omitted, \code{text} is selected for
    smaller data sets, and \code{points} for larger.}
  \item{xlim, ylim}{the x and y limits (min,max) of the plot.}
  \item{labels}{Optional text to be used instead of row names. If you
    use this, it is good to check the default labels and their order
    using \code{labels} command.}
  \item{arrow.mul}{Factor to expand arrows in the graph.  Arrows will be
    scaled automatically to fit the graph if this is missing.}
  \item{head.arrow}{Default length of arrow heads.}
  \item{select}{Items to be displayed.  This can either be a logical
    vector which is \code{TRUE} for displayed items or a vector of indices
    of displayed items.}
  \item{const}{General scaling constant to \code{rda} scores. The
    default is to use a constant that gives biplot scores, that is,
    scores that approximate original data (see \code{\link{vignette}}
    on \sQuote{Design Decisions} with \code{browseVignettes("vegan")}
    for details and discussion). If \code{const} is a vector of two
    items, the first is used for species, and the second item for site
    scores.}
  \item{droplist}{Return a matrix instead of a named list when only one
    kind of scores were requested.}
  \item{axis.bp}{Draw \code{\link{axis}} for biplot arrows.}
  \item{digits}{Number of digits in output.}

  \item{...}{Parameters passed to graphical functions. These will be
    applied to all score types, but will be superseded by score type
    parameters list (except \code{type = "none"} which will only draw
    the frame).}
}

\details{
  Same \code{plot} function will be used for \code{\link{cca}} and
  \code{\link{rda}}. This produces a quick, standard plot with current
  \code{scaling}.

  The \code{plot} function sets colours (\code{col}), plotting
  characters (\code{pch}) and character sizes (\code{cex}) to
  default values. The defaults can be changed with global parameters
  (\dQuote{dot arguments}) applied to all score types, or a list of
  parameters for a specified score type which take precedence over
  global parameters and defaults. If you use \code{points} or
  \code{text} to add scores to an existing plot, you must be careful to
  use exactly same \code{scaling}, \code{choices} and other parameters.

  The \code{plot} function returns (invisible) \code{\link{ordiplot}}
  object. You can save this object and use it to construct your plot
  with \code{ordiplot} functions \code{points} and \code{text}. These
  functions can be used in pipe (\code{|>}) which allows incremental
  building of plots with full control of graphical parameters for each
  score type. With pipe it is best to first create an empty plot with
  \code{plot(<cca-result>, type = "n")} and then add elements with
  \code{points}, \code{text} of \code{\link{ordiplot}} or
  \code{\link{ordilabel}}. Within pipe, the first argument should be a
  quoted score type, and then the grapcical parameters. The full
  object may contain scores with names
  \Sexpr[results=rd,stage=build]{require(vegan, quietly=TRUE);
    data(dune, dune.env);
    noquote(paste0(sQuote(names(scores(rda(dune ~ ., dune.env)))),
    collapse=", "))}
  (some of these may be missing depending on your model). The first
  \code{plot} will set the dimensions of graph, and if you do not use
  some score type there may be empty white space.

  Environmental variables receive a special treatment. With
  \code{display="bp"}, arrows will be drawn. These are labelled with
  \code{text} and unlabelled with \code{points}. The arrows have
  basically unit scaling, but if sites were scaled (\code{scaling}
  \code{"sites"} or \code{"symmetric"}), the scores of requested axes
  are adjusted to the plotting area.  With
  \code{scaling = "species"} or \code{scaling = "none"}, the arrows will
  be consistent with vectors fitted to linear combination scores
  (\code{display = "lc"} in function \code{\link{envfit}}), but with
  other scaling alternatives they will differ. The basic \code{plot}
  function uses a simple heuristics for adjusting the unit-length arrows
  to the current plot area, but the user can give the expansion factor
  in \code{arrow.mul}.  With \code{display="cn"} the centroids of levels
  of \code{\link{factor}} variables are displayed. With this option continuous
  variables still are presented as arrows and ordered factors as arrows
  and centroids. With \code{display = "reg"} arrows will be drawn for
  regression coefficients (a.k.a. canonical coefficients) of constraints
  and conditions. Biplot arrows can be interpreted individually, but
  regression coefficients must be interpreted all together: the LC score
  for each site is the sum of regressions displayed by arrows. The
  partialled out conditions are zero and not shown in biplot arrows, but
  they are shown for regressions, and show the effect that must be
  partialled out to get the LC scores. The biplot arrows are more
  standard and more easily interpreted, and regression arrows should be
  used only if you know that you need them.

  Palmer (1993) suggested using linear constraints (\dQuote{LC scores})
  in ordination diagrams, because these gave better results in
  simulations and site scores (\dQuote{WA scores}) are a step from
  constrained to unconstrained analysis.  However, McCune (1997) showed
  that noisy environmental variables (and all environmental measurements
  are noisy) destroy \dQuote{LC scores} whereas \dQuote{WA scores} were
  little affected.  Therefore the \code{plot} function uses site scores
  (\dQuote{WA scores}) as the default. This is consistent with the usage
  in statistics and other functions in \R (\code{\link[MASS]{lda}},
  \code{\link{cancor}}).  } \value{ The \code{plot} function returns
  invisibly a plotting structure which can be used by function
  \code{\link{identify.ordiplot}} to identify the points or other
  functions in the \code{\link{ordiplot}} family.  }

\author{Jari Oksanen }

\seealso{
  The function builds upon \code{\link{ordiplot}} and its \code{text}
  and \code{points} functions. See these to find available graphical
  parameters. Function \code{\link{ordilabel}} can also be used within a
  pipe.
}

\examples{
data(dune, dune.env)
mod <- cca(dune ~ A1 + Moisture + Management, dune.env)

## default and modified plot
plot(mod, scaling="sites")
plot(mod, scaling="sites", type = "text",
  sit.par = list(type = "points", pch=21, col="red", bg="yellow", cex=1.2),
  spe.par = list(col="blue", cex=0.8),
  cen.par = list(bg="white"))

## same with pipe
plot(mod, type="n", scaling="sites") |>
  points("sites", pch=21, col="red", bg = "yellow", cex=1.2) |>
  text("species", col="blue", cex=0.8) |>
  text("biplot") |>
  text("centroids", bg="white")

## catch the invisible result and use ordiplot support - the example
## will make a biplot with arrows for species and correlation scaling
pca <- rda(dune)
pl <- plot(pca, type="n", scaling="sites", correlation=TRUE)
with(dune.env, points(pl, "site", pch=21, col=1, bg=Management))
text(pl, "sp", arrow=TRUE, length=0.05, col=4, cex=0.6, xpd=TRUE)
with(dune.env, legend("bottomleft", levels(Management), pch=21,
   pt.bg=1:4, bty="n"))

## Pipe
plot(pca, type="n", scaling="sites", correlation=TRUE) |>
    points("sites", pch=21, col = 1, cex=1.5, bg = dune.env$Management) |>
    text("species", col = "blue", arrows = TRUE, xpd = TRUE, font = 3)

## Scaling can be numeric or more user-friendly names
## e.g. Hill's scaling for (C)CA
scrs <- scores(mod, scaling = "sites", hill = TRUE)
## or correlation-based scores in PCA/RDA
scrs <- scores(rda(dune ~ A1 + Moisture + Management, dune.env),
               scaling = "sites", correlation = TRUE)
}
\keyword{hplot}
\keyword{aplot}

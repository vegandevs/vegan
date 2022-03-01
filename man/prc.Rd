\name{prc}
\alias{prc}
\alias{summary.prc}
\alias{plot.prc}

\title{Principal Response Curves for Treatments with Repeated Observations }

\description{ Principal Response Curves (PRC) are a special case of
  Redundancy Analysis (\code{\link{rda}}) for multivariate responses in
  repeated observation design. They were originally suggested for
  ecological communities. They should be easier to interpret than
  traditional constrained ordination. They can also be used to study how
  the effects of a factor \code{A} depend on the levels of a factor
  \code{B}, that is \code{A + A:B}, in a multivariate response
  experiment.
}

\usage{
prc(response, treatment, time, ...)
\method{summary}{prc}(object, axis = 1, scaling = "symmetric", const,
        digits = 4, correlation = FALSE, ...)
\method{plot}{prc}(x, species = TRUE, select, scaling = "symmetric",
     axis = 1, correlation = FALSE, const, type = "l", xlab, ylab, ylim,
     lty = 1:5, col = 1:6, pch, legpos, cex = 0.8, ...)  
}

\arguments{
  \item{response}{Multivariate response data. Typically these are
    community (species) data. If the data are counts, they probably
    should be log transformed prior to the analysis. }
  \item{treatment}{A factor for treatments. }
  \item{time}{ An unordered factor defining the observations times in
    the repeated design.}
  \item{object, x}{An \code{prc} result object.}
  \item{axis}{Axis shown (only one axis can be selected).}
  \item{scaling}{Scaling of species scores, identical to the
    \code{scaling} in \code{\link{scores.rda}}.

    The type of scores can also be specified as one of \code{"none"},
    \code{"sites"}, \code{"species"}, or \code{"symmetric"}, which
    correspond to the values \code{0}, \code{1}, \code{2}, and \code{3}
    respectively. Argument \code{correlation} can be used in combination
    with these character descriptions to get the corresponding negative
    value.
  }
  \item{const}{General scaling constant for species scores (see
    \code{\link{scores.rda}} for details). Lower values will reduce the
    range of species scores, but will not influence the regression
    coefficients.}
  \item{digits}{Number of significant digits displayed.}
  \item{correlation}{logical; if \code{scaling} is a character
    description of the scaling type, \code{correlation} can be used to
    select correlation-like scores for PCA. See argument \code{scaling}
    for details.}
  \item{species}{Display species scores.}
  \item{select}{Vector to select displayed species. This can be a vector
    of indices or a logical vector which is \code{TRUE} for the selected
    species}
  \item{type}{Type of plot: \code{"l"} for lines, \code{"p"} for points
    or \code{"b"} for both.}
  \item{xlab, ylab}{Text to replace default axis labels.}
  \item{ylim}{Limits for the vertical axis.}
  \item{lty, col, pch}{Line type, colour and plotting characters
    (defaults supplied).}
  \item{legpos}{The position of the \code{\link{legend}}. A guess is
    made if this is not supplied, and \code{NA} will suppress legend. }
  \item{cex}{Character expansion for symbols and species labels.}
  \item{\dots}{ Other parameters passed to functions. }
}

\details{ PRC is a special case of \code{\link{rda}} with a single
  factor for \code{treatment} and a single factor for \code{time} points
  in repeated observations. In \pkg{vegan}, the corresponding
  \code{\link{rda}} model is defined as \code{rda(response ~ treatment *
  time + Condition(time))}. Since the \code{time} appears twice in the
  model formula, its main effects will be aliased, and only the main
  effect of treatment and interaction terms are available, and will be
  used in PRC. Instead of usual multivariate ordination diagrams, PRC
  uses canonical (regression) coefficients and species scores for a
  single axis. All that the current functions do is to provide a special
  \code{summary} and \code{plot} methods that display the
  \code{\link{rda}} results in the PRC fashion. The current version only
  works with default contrasts (\code{\link{contr.treatment}}) in which
  the coefficients are contrasts against the first level, and the levels
  must be arranged so that the first level is the control (or a
  baseline). If necessary, you must change the baseline level with
  function \code{\link{relevel}}.  
  

  Function \code{summary} prints the species scores and the
  coefficients. Function \code{plot} plots coefficients against
  \code{time} using \code{\link{matplot}}, and has similar defaults.
  The graph (and PRC) is meaningful only if the first \code{treatment}
  level is the control, as the results are contrasts to the first level
  when unordered factors are used. The plot also displays species scores
  on the right vertical axis using function
  \code{\link{linestack}}. Typically the number of species is so high
  that not all can be displayed with the default settings, but users can
  reduce character size or padding (\code{air}) in
  \code{\link{linestack}}, or \code{select} only a subset of the
  species. A legend will be displayed unless suppressed with
  \code{legpos = NA}, and the functions tries to guess where to put the
  legend if \code{legpos} is not supplied.
  
}
\value{
  The function is a special case of \code{\link{rda}} and returns its
  result object (see \code{\link{cca.object}}). However, a special
  \code{summary} and \code{plot} methods display returns differently
  than in \code{\link{rda}}.
}
\references{
  van den Brink, P.J. & ter Braak, C.J.F. (1999). Principal response
  curves: Analysis of time-dependent multivariate responses of
  biological community to stress. Environmental Toxicology and
  Chemistry, 18, 138--148.
}
\author{ Jari Oksanen and Cajo ter Braak}

\section{Warning }{The first level of \code{treatment} must be the
  control: use function \code{\link{relevel}} to guarantee the correct
  reference level. The current version will ignore user setting of
  \code{\link{contrasts}} and always use treatment contrasts
  (\code{\link{contr.treatment}}). The \code{time} must be an unordered
  factor.  }

\seealso{\code{\link{rda}}, \code{\link{anova.cca}}.}
\examples{
## Chlorpyrifos experiment and experimental design: Pesticide
## treatment in ditches (replicated) and followed over from 4 weeks
## before to 24 weeks after exposure 
data(pyrifos)
week <- gl(11, 12, labels=c(-4, -1, 0.1, 1, 2, 4, 8, 12, 15, 19, 24))
dose <- factor(rep(c(0.1, 0, 0, 0.9, 0, 44, 6, 0.1, 44, 0.9, 0, 6), 11))
ditch <- gl(12, 1, length=132)
# PRC
mod <- prc(pyrifos, dose, week)
mod            # RDA
summary(mod)   # PRC
logabu <- colSums(pyrifos)
plot(mod, select = logabu > 100)
## Ditches are randomized, we have a time series, and are only
## interested in the first axis
ctrl <- how(plots = Plots(strata = ditch,type = "free"),
    within = Within(type = "series"), nperm = 99)
anova(mod, permutations = ctrl, first=TRUE)
}
\keyword{ multivariate }


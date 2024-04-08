\name{ordiresids}
\alias{ordiresids}
 
\title{Plots of Residuals and Fitted Values for Constrained Ordination}

\description{
  The function provides \code{\link{plot.lm}} style diagnostic plots
  for the results of constrained ordination from \code{\link{cca}},
  \code{\link{rda}}, \code{\link{dbrda}} and
  \code{\link{capscale}}. Normally you do not need these plots,
  because ordination is descriptive and does not make assumptions on
  the distribution of the residuals. However, if you permute residuals
  in significance tests (\code{\link{anova.cca}}), you may be
  interested in inspecting that the residuals really are exchangeable
  and independent of fitted values.
}

\usage{
ordiresids(x, kind = c("residuals", "scale", "qqmath"),
   residuals = "working", type = c("p", "smooth", "g"),
   formula, ...)
}
\arguments{
  \item{x}{Ordination result from \code{\link{cca}}, \code{\link{rda}},
    \code{\link{dbrda}}, \code{\link{capscale}}. }
  \item{kind}{The type of plot: \code{"residuals"} plot residuals
    against fitted values, \code{"scale"} the square root of absolute
    residuals against fitted values, and \code{"qqmath"} the residuals
    against expected distribution (defaults \code{\link{qnorm}}),
    unless defined differently in the \code{formula}  argument. }
  \item{residuals}{The kind of residuals and fitted values, with alternatives
    \code{"working"}, \code{"response"}, \code{"standardized"} and
    \code{"studentized"} (see Details).}
  \item{type}{The type of plot. The argument is passed on to
    \pkg{lattice} functions. }
  \item{formula}{Formula to override the default plot. The formula can
    contain items \code{Fitted}, \code{Residuals}, \code{Species} and
    \code{Sites} (provided that names of species and sites are available
    in the ordination result).  }
  \item{\dots}{Other arguments passed to \pkg{lattice} functions. }
}

\details{
  The default plots are similar as in \code{\link{plot.lm}}, but they
  use \code{\link[lattice]{Lattice}} functions
  \code{\link[lattice]{xyplot}} and \code{\link[lattice]{qqmath}}. The
  alternatives have default formulae but these can be replaced by the
  user. The elements available in formula or in the \code{groups} argument
  are \code{Fitted}, \code{Residuals}, \code{Species} and \code{Sites}.

  With \code{residuals = "response"} and \code{residuals = "working"}
  the fitted values and residuals are found with functions
  \code{\link{fitted.cca}} and \code{\link{residuals.cca}}. With
  \code{residuals = "standardized"} the residuals are found with
  \code{\link{rstandard.cca}}, and with \code{residuals = "studentized"}
  they are found with \code{\link{rstudent.cca}}, and in both cases the
  fitted values are standardized with \code{\link{sigma.cca}}.
}

\value{
  The function returns a \code{\link[lattice]{Lattice}} object that can
  displayed as plot.
}

\author{Jari Oksanen }
\seealso{

  \code{\link{plot.lm}}, \code{\link{fitted.cca}},
  \code{\link{residuals.cca}}, \code{\link{rstandard.cca}},
  \code{\link{rstudent.cca}}, \code{\link{sigma.cca}},
  \code{\link[lattice]{Lattice}}, \code{\link[lattice]{xyplot}},
  \code{\link[lattice]{qqmath}}.  }


\examples{
data(varespec)
data(varechem)
mod <- cca(varespec ~ Al + P + K, varechem)
ordiresids(mod)
ordiresids(mod, formula = Residuals ~ Fitted | Species, residuals="standard",
   cex = 0.5)
}
\keyword{ hplot }


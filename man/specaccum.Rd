\name{specaccum}
\alias{specaccum}
\alias{print.specaccum}
\alias{summary.specaccum}
\alias{plot.specaccum}
\alias{lines.specaccum}
\alias{boxplot.specaccum}
\alias{fitspecaccum}
\alias{plot.fitspecaccum}
\alias{lines.fitspecaccum}
\alias{predict.specaccum}
\alias{predict.fitspecaccum}
\alias{AIC.fitspecaccum}
\alias{deviance.fitspecaccum}
\alias{specslope}

\title{Species Accumulation Curves}
\description{
  Function \code{specaccum} finds species accumulation curves or the
  number of species for a certain number of sampled sites or
  individuals. 
}
\usage{
specaccum(comm, method = "exact", permutations = 100,
          conditioned =TRUE, gamma = "jack1",  w = NULL, subset, ...)
\method{plot}{specaccum}(x, add = FALSE, random = FALSE, ci = 2, 
    ci.type = c("bar", "line", "polygon"), col = par("fg"), lty = 1,
    ci.col = col, ci.lty = 1, xlab, ylab = x$method, ylim,
    xvar = c("sites", "individuals", "effort"), ...)
\method{boxplot}{specaccum}(x, add = FALSE, ...)
fitspecaccum(object, model, method = "random", ...)
\method{plot}{fitspecaccum}(x, col = par("fg"), lty = 1, xlab = "Sites", 
    ylab = x$method, ...) 
\method{predict}{specaccum}(object, newdata, interpolation = c("linear", "spline"), ...)
\method{predict}{fitspecaccum}(object, newdata, ...)
specslope(object, at)
}

\arguments{
  \item{comm}{Community data set.}
  \item{method}{Species accumulation method (partial match). Method
    \code{"collector"}
    adds sites in the order they happen to be in the data,
    \code{"random"} adds sites in random order, \code{"exact"} finds the
    expected (mean) species richness, \code{"coleman"} finds the
    expected richness following
    Coleman et al. 1982, and \code{"rarefaction"} finds the mean when
    accumulating individuals instead of sites.  }
  \item{permutations}{Number of permutations with \code{method = "random"}.
    Usually an integer giving the number permutations, but can also be a
    list of control values for the permutations as returned by the
    function \code{\link[permute]{how}}, or a permutation matrix where
    each row gives the permuted indices.
  }
  \item{conditioned}{ Estimation of standard deviation is conditional on
    the empirical dataset for the exact SAC}
  \item{gamma}{Method for estimating the total extrapolated number of species in the
    survey area by function \code{\link{specpool}}}
  \item{w}{Weights giving the sampling effort.}
  \item{subset}{logical expression indicating sites (rows) to keep: missing
    values are taken as \code{FALSE}.}
  \item{x}{A \code{specaccum} result object}
  \item{add}{Add to an existing graph.}
  \item{random}{Draw each random simulation separately instead of
    drawing their average and confidence intervals.}
  \item{ci}{Multiplier used to get confidence intervals from standard
    deviation (standard error of the estimate). Value \code{ci = 0}
    suppresses drawing confidence intervals.}
  \item{ci.type}{Type of confidence intervals in the graph: \code{"bar"}
    draws vertical bars, \code{"line"} draws lines, and
    \code{"polygon"} draws a shaded area.}
  \item{col}{Colour for drawing lines.}
  \item{lty}{line type (see \code{\link{par}}).}
  \item{ci.col}{Colour for drawing lines or filling the
    \code{"polygon"}.}
  \item{ci.lty}{Line type for confidence intervals or border of the
    \code{"polygon"}.}
  \item{xlab,ylab}{Labels for \code{x} (defaults \code{xvar}) and
    \code{y} axis.}
  \item{ylim}{the y limits of the plot.}
  \item{xvar}{Variable used for the horizontal axis:
    \code{"individuals"} can be used only with
    \code{method = "rarefaction"}. }
  \item{object}{Either a community data set or fitted \code{specaccum} model.}
  \item{model}{Nonlinear regression model (\code{\link{nls}}). See Details.}
  
  \item{newdata}{Optional data used in prediction interpreted as
    number of sampling units (sites). If missing, fitted values are
    returned. }
  
  \item{interpolation}{Interpolation method used with \code{newdata}.}

  \item{at}{Number of plots where the slope is evaluated. Can be a
    real number.}

  \item{...}{Other parameters to functions.}
}
\details{

  Species accumulation curves (SAC) are used to compare diversity
  properties of community data sets using different accumulator
  functions. The classic method is \code{"random"} which finds the mean
  SAC and its standard deviation from random permutations of the data,
  or subsampling without replacement (Gotelli & Colwell 2001).  The
  \code{"exact"} method finds the expected SAC using sample-based
  rarefaction method that has been independently developed numerous
  times (Chiarucci et al. 2008) and it is often known as Mao Tau
  estimate (Colwell et al. 2012).  The unconditional standard deviation
  for the exact SAC represents a moment-based estimation that is not
  conditioned on the empirical data set (sd for all samples > 0). The
  unconditional standard deviation is based on an estimation of the
  extrapolated number of species in the survey area (a.k.a. gamma
  diversity), as estimated by function \code{\link{specpool}}. The
  conditional standard deviation that was developed by Jari Oksanen (not
  published, sd=0 for all samples). Method \code{"coleman"} finds the
  expected SAC and its standard deviation following Coleman et
  al. (1982).  All these methods are based on sampling sites without
  replacement. In contrast, the \code{method = "rarefaction"} finds the
  expected species richness and its standard deviation by sampling
  individuals instead of sites.  It achieves this by applying function
  \code{\link{rarefy}} with number of individuals corresponding to
  average number of individuals per site.

  Methods \code{"random"} and \code{"collector"} can take weights
  (\code{w}) that give the sampling effort for each site.  The weights
  \code{w} do not influence the order the sites are accumulated, but
  only the value of the sampling effort so that not all sites are
  equal. The summary results are expressed against sites even when the
  accumulation uses weights (methods \code{"random"},
  \code{"collector"}), or is based on individuals
  (\code{"rarefaction"}).  The actual sampling effort is given as item
  \code{Effort} or \code{Individuals} in the printed result. For
  weighted \code{"random"} method the effort refers to the average
  effort per site, or sum of weights per number of sites. With
  weighted \code{method = "random"}, the averaged species richness is
  found from linear interpolation of single random permutations.
  Therefore at least the first value (and often several first) have
  \code{NA} richness, because these values cannot be interpolated in
  all cases but should be extrapolated.  The \code{plot} function
  defaults to display the results as scaled to sites, but this can be
  changed selecting \code{xvar = "effort"} (weighted methods) or
  \code{xvar = "individuals"} (with \code{method = "rarefaction"}).
 
  The \code{summary} and \code{boxplot} methods are available for
  \code{method = "random"}.

  Function \code{predict} for \code{specaccum} can return the values
  corresponding to \code{newdata}. With \code{method} \code{"exact"},
  \code{"rarefaction"} and \code{"coleman"} the function uses analytic
  equations for interpolated non-integer values, and for other methods
  linear (\code{\link{approx}}) or spline (\code{\link{spline}})
  interpolation. If \code{newdata} is not given, the function returns
  the values corresponding to the data. NB., the fitted values with
  \code{method="rarefaction"} are based on rounded integer counts, but
  \code{predict} can use fractional non-integer counts with
  \code{newdata} and give slightly different results.

  Function \code{fitspecaccum} fits a nonlinear (\code{\link{nls}})
  self-starting species accumulation model. The input \code{object}
  can be a result of \code{specaccum} or a community in data frame. In
  the latter case the function first fits a \code{specaccum} model and
  then proceeds with fitting the nonlinear model. The function can
  apply a limited set of nonlinear regression models suggested for
  species-area relationship (Dengler 2009). All these are
  \code{\link{selfStart}} models. The permissible alternatives are
  \code{"arrhenius"} (\code{\link{SSarrhenius}}), \code{"gleason"}
  (\code{\link{SSgleason}}), \code{"gitay"} (\code{\link{SSgitay}}),
  \code{"lomolino"} (\code{\link{SSlomolino}}) of \pkg{vegan}
  package. In addition the following standard \R models are available:
  \code{"asymp"} (\code{\link{SSasymp}}), \code{"gompertz"}
  (\code{\link{SSgompertz}}), \code{"michaelis-menten"})
  (\code{\link{SSmicmen}}), \code{"logis"} (\code{\link{SSlogis}}),
  \code{"weibull"} (\code{\link{SSweibull}}). See these functions for
  model specification and details. 

  When weights \code{w} were used the fit is based on accumulated
  effort and in \code{model = "rarefaction"} on accumulated number of
  individuals.  The \code{plot} is still based on sites, unless other
  alternative is selected with \code{xvar}.

  Function \code{predict} for \code{fitspecaccum} uses
  \code{\link{predict.nls}}, and you can pass all arguments to that
  function. In addition, \code{fitted}, \code{residuals} and
  \code{coef}, \code{AIC} and \code{deviance} work on the result
  object.

  Function \code{specslope} evaluates the derivative of the species
  accumulation curve at given number of sample plots, and gives the
  rate of increase in the number of species. The function works with
  \code{specaccum} result object when this is based on analytic models
  \code{"exact"}, \code{"rarefaction"} or \code{"coleman"}, and with
  non-linear regression results of \code{fitspecaccum}.

  Nonlinear regression may fail for any reason, and some of the
  \code{fitspecaccum} models are fragile and may not succeed.  
}

\value{ Function \code{specaccum} returns an object of class
  \code{"specaccum"}, and \code{fitspecaccum} a model of class
  \code{"fitspecaccum"} that adds a few items to the
  \code{"specaccum"} (see the end of the list below):

  \item{call }{Function call.}
  \item{method}{Accumulator method.}
  \item{sites}{Number of sites.  For \code{method = "rarefaction"} this
    is the number of sites corresponding to a certain number of
    individuals and generally not an integer, and the average
    number of individuals is also returned in item \code{individuals}.} 
  \item{effort}{Average sum of weights corresponding to the number of
    sites when model was fitted with argument \code{w}}
  \item{richness}{The number of species corresponding to number of
    sites.  With \code{method = "collector"} this is the observed
    richness, for other methods the average or expected richness.}
  \item{sd}{The standard deviation of SAC (or its standard error). This
    is \code{NULL} in \code{method = "collector"}, and it
    is estimated from permutations in \code{method = "random"}, and from
    analytic equations in other methods.}
  \item{perm}{Permutation results with \code{method = "random"} and
    \code{NULL} in other cases. Each column in \code{perm} holds one
    permutation.}
  \item{weights}{Matrix of accumulated weights corresponding to the
    columns of the \code{perm} matrix when model was fitted with
    argument \code{w}.}
  \item{fitted, residuals, coefficients}{Only in \code{fitspecacum}:
     fitted values, residuals and nonlinear model coefficients. For
     \code{method = "random"} these are matrices with a column for
     each random accumulation.}

  \item{models}{Only in \code{fitspecaccum}: list of fitted
    \code{\link{nls}} models (see Examples on accessing these models).} 

}

\references{
  Chiarucci, A., Bacaro, G., Rocchini, D. & Fattorini,
  L. (2008). Discovering and rediscovering the sample-based rarefaction
  formula in the ecological literature. \emph{Commun. Ecol.} 9:
  121--123.

  Coleman, B.D, Mares, M.A., Willis, M.R. & Hsieh,
  Y. (1982). Randomness, area and species richness. \emph{Ecology} 63:
  1121--1133. 

  Colwell, R.K., Chao, A., Gotelli, N.J., Lin, S.Y., Mao, C.X., Chazdon,
  R.L. & Longino, J.T. (2012). Models and estimators linking
  individual-based and sample-based rarefaction, extrapolation and
  comparison of assemblages. \emph{J. Plant Ecol.} 5: 3--21.

  Dengler, J. (2009). Which function describes the species-area
  relationship best? A review and empirical evaluation. 
  \emph{Journal of Biogeography} 36, 728--744.

  Gotelli, N.J. & Colwell, R.K. (2001). Quantifying biodiversity:
  procedures and pitfalls in measurement and comparison of species
  richness. \emph{Ecol. Lett.} 4, 379--391.
}
\author{Roeland Kindt \email{r.kindt@cgiar.org} and Jari Oksanen.}
\note{
  The SAC with \code{method = "exact"} was
  developed by Roeland Kindt, and its standard deviation by Jari
  Oksanen (both are unpublished). The \code{method = "coleman"}
  underestimates the SAC because it does not handle properly sampling
  without replacement.  Further, its standard deviation does not take
  into account species correlations, and is generally too low. }

\seealso{\code{\link{rarefy}} and \code{\link{rrarefy}} are related
  individual based models. Other accumulation models are
  \code{\link{poolaccum}} for extrapolated richness, and
  \code{\link{renyiaccum}} and \code{\link{tsallisaccum}} for
  diversity indices.  Underlying graphical functions are
  \code{\link{boxplot}}, \code{\link{matlines}},
  \code{\link{segments}} and \code{\link{polygon}}.
}
\examples{
data(BCI)
sp1 <- specaccum(BCI)
sp2 <- specaccum(BCI, "random")
sp2
summary(sp2)
plot(sp1, ci.type="poly", col="blue", lwd=2, ci.lty=0, ci.col="lightblue")
boxplot(sp2, col="yellow", add=TRUE, pch="+")
## Fit Lomolino model to the exact accumulation
mod1 <- fitspecaccum(sp1, "lomolino")
coef(mod1)
fitted(mod1)
plot(sp1)
## Add Lomolino model using argument 'add'
plot(mod1, add = TRUE, col=2, lwd=2)
## Fit Arrhenius models to all random accumulations
mods <- fitspecaccum(sp2, "arrh")
plot(mods, col="hotpink")
boxplot(sp2, col = "yellow", border = "blue", lty=1, cex=0.3, add= TRUE)
## Use nls() methods to the list of models
sapply(mods$models, AIC)
}
\keyword{univar}
\keyword{models}


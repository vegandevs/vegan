\name{specaccum}
\alias{specaccum}
\alias{print.specaccum}
\alias{summary.specaccum}
\alias{plot.specaccum}
\alias{boxplot.specaccum}

\title{Species Accumulation Curves }
\description{
  Function \code{specaccum} finds species accumulation curves or the
  number of species for a certain number of sampled sites or
  individuals. 
}
\usage{
specaccum(comm, method = "exact", permutations = 100,
          conditioned =TRUE, gamma = "jack1",  ...)
\method{plot}{specaccum}(x, add = FALSE, ci = 2, ci.type = c("bar", "line", "polygon"), 
    col = par("fg"), ci.col = col, ci.lty = 1, xlab = "Sites", 
    ylab = x$method, ylim, ...)
\method{boxplot}{specaccum}(x, add = FALSE, ...)
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
  \item{permutations}{Number of permutations with \code{method =
      "random"}.}
  \item{conditioned}{ Estimation of standard deviation is conditional on
    the empirical dataset for the exact SAC}
  \item{gamma}{Method for estimating the total extrapolated number of species in the
    survey area by function \code{\link{specpool}}} 
  \item{x}{A \code{specaccum} result object}
  \item{add}{Add to an existing graph.}
  \item{ci}{Multiplier used to get confidence intervals from standard
    deviation (standard error of the estimate). Value \code{ci = 0}
    suppresses drawing confidence intervals.}
  \item{ci.type}{Type of confidence intervals in the graph: \code{"bar"}
    draws vertical bars, \code{"line"} draws lines, and
    \code{"polygon"} draws a shaded area.}
  \item{col}{Colour for drawing lines.}
  \item{ci.col}{Colour for drawing lines or filling the
    \code{"polygon"}.}
  \item{ci.lty}{Line type for confidence intervals or border of the
    \code{"polygon"}.}
  \item{xlab,ylab}{Labels for \code{x} and \code{y} axis.}
  \item{ylim}{the y limits of the plot.}
  \item{...}{Other parameters to functions.}
}
\details{
  Species accumulation curves (SAC) are used to compare diversity properties
  of community data sets using different accumulator functions. The
  classic method is \code{"random"} which finds the mean SAC and its
  standard deviation from random permutations of the data, or
  subsampling without replacement (Gotelli & Colwell 2001).
  The \code{"exact"} method finds the
  expected SAC using the method that was independently developed by
  Ugland et al. (2003), Colwell et al. (2004) and Kindt et al. (2006). 
  The unconditional standard deviation for the exact SAC represents a
  moment-based estimation that is not conditioned on the empirical data
  set (sd for all samples > 0), unlike the conditional standard deviation
  that was developed by Jari Oksanen (not published, sd=0 for all
  samples). The unconditional standard deviation is based on an estimation
  of the total extrapolated number of species in the survey area
  (a.k.a. gamma diversity), as estimated by
  function \code{\link{specpool}}.
  Method \code{"coleman"} finds the expected SAC and its standard
  deviation following Coleman et al. (1982).  All these methods are
  based on sampling sites without replacement. In contrast, the
  \code{method = "rarefaction"} finds the expected species richness and
  its standard deviation by sampling individuals instead of sites. It
  achieves this by applying function \code{\link{rarefy}} with number of individuals
  corresponding to average number of individuals per site.

  The function has a \code{plot} method. In addition, \code{method =
    "random"} has \code{summary} and \code{boxplot} methods. 
}

\value{
  The function returns an object of class \code{"specaccum"} with items:
  \item{call }{Function call.}
  \item{method}{Accumulator method.}
  \item{sites}{Number of sites.  For \code{method = "rarefaction"} this
    is the number of sites corresponding to a certain number of
    individuals and generally not an integer, and the average
    number of individuals is also retunred in item \code{individuals}.} 
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
}
\references{
  Coleman, B.D, Mares, M.A., Willis, M.R. & Hsieh,
  Y. (1982). Randomness, area and species richness. \emph{Ecology} 63:
  1121--1133. 

  Colwell, R.K., Mao, C.X. & Chang, J. (2004). Interpolating,
  extrapolating, and comparing incidence-based species accumulation
  curves. \emph{Ecology} 85: 2717--2727.  

  Gotellli, N.J. & Colwell, R.K. (2001). Quantifying biodiversity:
  procedures and pitfalls in measurement and comparison of species
  richness. \emph{Ecol. Lett.} 4, 379--391.

  Kindt, R. (2003). Exact species richness for sample-based accumulation
  curves. \emph{Manuscript.}
  
  Kindt R., Van Damme, P. & Simons, A.J. (2006) Patterns of species
  richness at varying scales in western Kenya: planning for
  agroecosystem diversification. \emph{Biodiversity and Conservation}, online
  first: DOI 10.1007/s10531-005-0311-9 

  Ugland, K.I., Gray, J.S. & Ellingsen, K.E. (2003). The
  species-accumulation curve and estimation of species richness. \emph{Journal
  of Animal Ecology} 72: 888--897.
  
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
  \code{\link{poolaccum}} for extrapoltated richenss, and
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
}
\keyword{univar }


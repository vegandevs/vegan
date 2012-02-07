\encoding{UTF-8}
\name{renyi}
\alias{renyi}
\alias{plot.renyi}
\alias{renyiaccum}
\alias{plot.renyiaccum}
\alias{persp.renyiaccum}
\alias{rgl.renyiaccum}

\title{Renyi and Hill Diversities and Corresponding Accumulation Curves }
\description{
  Function \code{renyi}  find \enc{Rényi}{Renyi} diversities with any
  scale or the corresponding Hill number (Hill 1973).  Function
  \code{renyiaccum} finds these statistics with accumulating sites.
}
\usage{
renyi(x, scales = c(0, 0.25, 0.5, 1, 2, 4, 8, 16, 32, 64, Inf), hill = FALSE)
\method{plot}{renyi}(x, ...)
renyiaccum(x, scales = c(0, 0.5, 1, 2, 4, Inf), permutations = 100, 
    raw = FALSE, ...)
\method{plot}{renyiaccum} (x, what = c("mean", "Qnt 0.025", "Qnt 0.975"), type = "l", 
    ...)
\method{persp}{renyiaccum} (x, theta = 220, col = heat.colors(100), zlim, ...)
rgl.renyiaccum(x, rgl.height = 0.2, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{Community data matrix or plotting object. }
  \item{scales}{Scales of \enc{Rényi}{Renyi} diversity.}
  \item{hill}{Calculate Hill numbers.}
  \item{permutations}{Number of random permutations in accumulating
    sites.}
  \item{raw}{if \code{FALSE} then return summary statistics of
    permutations, and if \code{TRUE} then returns the individual
    permutations.}
  \item{what}{Items to be plotted.}
  \item{type}{Type of plot, where \code{type = "l"} means lines.}
  \item{theta}{Angle defining the viewing direction (azimuthal) in
    \code{\link{persp}}.}
  \item{col}{Colours used for surface. Single colour will be passed on,
    and vector colours will be
    selected by the midpoint of a rectangle in \code{\link{persp}}. }
  \item{zlim}{Limits of vertical axis.}
  \item{rgl.height}{Scaling of vertical axis.}
  \item{\dots}{Other arguments which are passed to \code{renyi} and
    to graphical functions.}

}
\details{
  Common \code{\link{diversity}} indices are special cases of
  \enc{Rényi}{Renyi} diversity
  \deqn{H_a = \frac{1}{1-a} \log \sum p_i^a}{H.a = 1/(1-a) log
    sum(p^a)}
  where \eqn{a} is a scale parameter, and Hill (1975) suggested to
  use so-called ``Hill numbers'' defined as \eqn{N_a = \exp(H_a)}{N.a =
    exp(H.a)}.  Some Hill numbers are the number of species with
  \eqn{a = 0}, \eqn{\exp(H')}{exp(H')} or the exponent of Shannon
  diversity with \eqn{a = 1}, inverse Simpson with \eqn{a = 2} and
  \eqn{1/ \max(p_i)}{1/max(p)} with \eqn{a = \infty}{a = Inf}. According
  to the theory of diversity ordering, one community can be regarded as
  more diverse than another only if its \enc{Rényi}{Renyi} diversities are all higher 
  (\enc{Tóthmérész}{Tothmeresz}  1995).

  The \code{plot} method for \code{renyi} uses \pkg{lattice} graphics,
  and displays the diversity values against each scale in separate panel
  for each site together with minimum, maximum and median values in the
  complete data.

  Function \code{renyiaccum} is similar to \code{\link{specaccum}} but
  finds \enc{Rényi}{Renyi} or Hill diversities at given \code{scales}
  for random permutations of accumulated sites.  Its \code{plot}
  function uses \pkg{lattice} function \code{\link[lattice]{xyplot}} to
  display the accumulation curves for each value of \code{scales} in a
  separate panel.  In addition, it has a \code{persp} method to plot the
  diversity surface against scale and number and sites. Dynamic graphics
  with \code{rgl.renyiaccum} use \pkg{rgl} package, and produces similar
  surface as \code{persp} with a mesh showing the empirical confidence levels.
}
\value{
  Function \code{renyi} returns a data frame of selected
  indices.  Function \code{renyiaccum} with argument \code{raw = FALSE}
  returns a three-dimensional array, where the first dimension are the
  accumulated sites, second dimension are the diversity scales, and
  third dimension are the summary statistics \code{mean}, \code{stdev},
  \code{min}, \code{max}, \code{Qnt 0.025} and \code{Qnt 0.975}.  With
  argument \code{raw = TRUE} the statistics on the third dimension are
  replaced with individual permutation results. 
}
\references{
  \url{http://www.worldagroforestry.org/resources/databases/tree-diversity-analysis}

  Hill, M.O. (1973). Diversity and evenness: a unifying notation and its
   consequences. \emph{Ecology} 54, 427--473.

   Kindt R, Van Damme P, Simons AJ. 2006. Tree diversity in western
   Kenya: using profiles to characterise richness and 
   evenness. \emph{Biodiversity and Conservation} 15: 1253-1270.
   
   \enc{Tóthmérész}{Tothmeresz}, B. (1995). Comparison of different methods for diversity
   ordering. \emph{Journal of Vegetation Science} 6, 283--290.
}
\author{ Roeland Kindt \email{r.kindt@cgiar.org} and Jari Oksanen }

\seealso{\code{\link{diversity}} for diversity indices, and
  \code{\link{specaccum}} for ordinary species accumulation curves, and
  \code{\link[lattice]{xyplot}}, \code{\link{persp}} and
  \code{\link[rgl]{rgl}} for controlling graphics. } 
\examples{
data(BCI)
i <- sample(nrow(BCI), 12)
mod <- renyi(BCI[i,])
plot(mod)
mod <- renyiaccum(BCI[i,])
plot(mod, as.table=TRUE, col = c(1, 2, 2))
persp(mod)
} 
\keyword{univar }


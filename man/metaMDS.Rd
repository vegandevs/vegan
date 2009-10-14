\name{metaMDS}
\alias{metaMDS}
\alias{metaMDSdist}
\alias{metaMDSiter}
\alias{metaMDSredist}
\alias{initMDS}
\alias{postMDS}
\alias{print.metaMDS}
\alias{plot.metaMDS}
\alias{points.metaMDS}
\alias{text.metaMDS}
\alias{scores.metaMDS}
\alias{metaMDSrotate}

\title{Nonmetric Multidimensional Scaling with Stable Solution from
  Random Starts, Axis Scaling and Species Scores}
\description{
  Function \code{metaMDS} uses \code{\link[MASS]{isoMDS}} to perform
  Nonmetric Multidimensional Scaling (NMDS), but tries to find a stable
  solution using several random starts (function \code{initMDS}). In
  addition, it standardizes the scaling in
  the result, so that the configurations are easier
  to interpret (function \code{postMDS}), and adds species scores to the
  site ordination (function \code{\link{wascores}}). 
}
\usage{
metaMDS(comm, distance = "bray", k = 2, trymax = 20, autotransform =TRUE,
        noshare = 0.1, wascores = TRUE, expand = TRUE, trace = 1,
        plot = FALSE, previous.best, old.wa = FALSE, ...)
\method{plot}{metaMDS}(x, display = c("sites", "species"), choices = c(1, 2),
     type = "p", shrink = FALSE,  ...)
\method{points}{metaMDS}(x, display = c("sites", "species"),
        choices = c(1,2), shrink = FALSE, select, ...)
\method{text}{metaMDS}(x, display = c("sites", "species"), labels, 
        choices = c(1,2), shrink = FALSE, select, ...)
\method{scores}{metaMDS}(x, display = c("sites", "species"), shrink = FALSE, 
        choices, ...)
metaMDSdist(comm, distance = "bray", autotransform = TRUE, noshare = 0.1, 
    trace = 1, commname, zerodist = "fail", distfun = vegdist, ...)
metaMDSiter(dist, k = 2, trymax = 20, trace = 1, plot = FALSE, previous.best, 
    ...)   
initMDS(x, k=2)
postMDS(X, dist, pc=TRUE, center=TRUE, halfchange, threshold=0.8,
        nthreshold=10, plot=FALSE, ...)
metaMDSredist(object, ...)
metaMDSrotate(object, vec, choices, ...)
}

\arguments{
  \item{comm}{Community data. Alternatively, dissimilarities either as
    a \code{\link{dist}} structure or as a symmetric square matrix. 
    In the latter case all other stages are skipped except random 
    starts and centring and pc rotation of axes. }
  \item{distance}{Dissimilarity index used in \code{\link{vegdist}}.}
  \item{k}{Number of dimensions in \code{\link[MASS]{isoMDS}}.}
  \item{trymax}{Maximum number of random starts in search of stable
    solution.}
  \item{autotransform}{Use simple heuristics for possible data
    transformation (see below).}
  \item{noshare}{Proportion of site pairs with no shared species to
    trigger \code{\link{stepacross}} to find flexible shortest paths
    among dissimilarities.}
  \item{wascores}{Calculate species scores using function
    \code{\link{wascores}}.}
  \item{expand}{Expand weighted averages of species in
    \code{\link{wascores}}.}
  \item{trace}{Trace the function; \code{trace = 2} or higher will be
    more voluminous.}
  \item{plot}{Graphical tracing: plot interim results. You may want to set
    \code{par(ask = TRUE)} with this option.}
  \item{previous.best}{Start searches from a previous
    solutions. Otherwise use \code{\link[MASS]{isoMDS}} default for the
    starting solution.}
  \item{old.wa}{Use the old way of calculating WA scores for species:
    in vegan versions 1.12-5 and 1.11-2 WA scores were based on
    untransformed data even when data were transformed in analysis, but
    since then the similar transformation will be used in WA scores as
    in ordination.}
  \item{x}{Dissimilarity matrix for \code{\link[MASS]{isoMDS}} or \code{plot}
    object.} 
  \item{choices}{Axes shown.}
  \item{type}{Plot type: \code{"p"} for points, \code{"t"} for text, and
   \code{"n"} for axes only.}
 \item{display}{Display \code{"sites"} or \code{"species"}.}
 \item{shrink}{Shrink back species scores if they were expanded
   originally.}
 \item{labels}{Optional test to be used instead of row names.}
 \item{select}{Items to be displayed.  This can either be a logical
   vector which is \code{TRUE} for displayed items or a vector of indices
   of displayed items.}
 \item{X}{Configuration from multidimensional scaling. }
 \item{commname}{The name of \code{comm}: should not be given if the
   function is called directly.}
 \item{zerodist}{Handling of zero dissimilarities: either \code{"fail"}
   or \code{"add"} a small positive value, or \code{"ignore"}.}
 \item{distfun}{Dissimilarity function. Any function returning a
   \code{dist} object and accepting argument \code{method} can be used
   (but some extra arguments may cause name conflicts).}
 \item{dist}{Dissimilarity matrix used in multidimensional scaling. }
  \item{pc}{Rotate to principal components. }
  \item{center}{Centre the configuration. }
 \item{halfchange}{Scale axes to half-change units. This defaults
  \code{TRUE} when dissimilarities were evaluated within
  \code{metaMDS} and the dissimilarity index has an upper limit of
  \eqn{1}. If \code{FALSE}, the ordination dissimilarities are scaled
  to the same range as the input dissimilarities.}
 \item{threshold}{Largest dissimilarity used in half-change scaling. }
 \item{nthreshold}{ Minimum number of points in half-change scaling. }
 \item{object}{A result object from \code{metaMDS}.}
 \item{vec}{A continuous site variable (vector).}
 \item{...}{Other parameters passed to functions.}
}
\details{
  Non-metric Multidimensional Scaling (NMDS) is commonly regarded as the
  most robust unconstrained ordination method in community ecology (Minchin
  1987).  Functions \code{initMDS} and \code{postMDS} together with some
  other functions are intended to 
  help run NMDS wit \code{\link[MASS]{isoMDS}} like recommended by
  Minchin (1987).  Function \code{metaMDS} combines all recommendations
  into one command for a shotgun style analysis. The complete steps in
  \code{metaMDS} are:
  \enumerate{
    \item Transformation: If the data values are larger than common
    class scales, the function performs a Wisconsin double
    standardization using \code{\link{wisconsin}}.  If the values look
    very large, the function also performs \code{\link{sqrt}}
    transformation. Both of these standardization are generally found
    to improve the results. However, the limits are completely
    arbitrary (at present, data maximum 50 triggers \code{\link{sqrt}}
    and >9 triggers \code{\link{wisconsin}}). If you want to have a full
    control of the analysis, you
    should set \code{autotransform = FALSE} and make explicit
    standardization in the command.

    \item Choice of dissimilarity: For a good result, you should use
    dissimilarity indices that have a good rank order relation to
    ordering sites along gradients (Faith et al. 1987).
    The default is Bray dissimilarity, because it often is the test
    winner. However, any other
    dissimilarity index in \code{\link{vegdist}} can be used. Function
    \code{\link{rankindex}} can be used for finding the test winner for
    you data and gradients.

    \item Step-across dissimilarities: Ordination may be very difficult
    if a large proportion of sites have no shared species. In this case,
    the results may be improved with \code{\link{stepacross}}
    dissimilarities, or flexible shortest paths among all sites. The
    \code{\link{stepacross}} is triggered by option \code{noshare}. If
    you do not like manipulation of original distances, you should set
    \code{noshare = 1}.  

    \item NMDS with random starts: NMDS easily gets trapped into local
    optima, and you must start NMDS several times from random start to
    be confident that you have found the global solution. The default in
    \code{\link[MASS]{isoMDS}} is to start from
    metric scaling (with \code{\link{cmdscale}}) which typically is
    close to a local optimum. The strategy in \code{metaMDS} is to first
    run a default \code{\link[MASS]{isoMDS}}, or use the
    \code{previous.best} solution if supplied, and take
    its solution as the standard (\code{Run 0}). Then \code{metaMDS}
    starts \code{\link[MASS]{isoMDS}} from several
    random starts (maximum number is given by \code{trymax}). If a
    solution is better (has a lower stress) than the previous standard,
    it is taken as the new standard. If the solution is better or close
    to a standard, \code{metaMDS} compares two solutions using
    Procrustes analysis using function \code{\link{procrustes}} with
    option \code{symmetric = TRUE}. If the two solutions are very
    similar in their Procrustes \code{rmse} and the largest residual is
    very small, the solutions are regarded as convergent and the best
    one is saved. Please note that the conditions are stringent, and you
    may have found good and relatively stable solutions although the
    function is not yet satisfied. Setting \code{trace = TRUE} will
    monitor the final stresses, and \code{plot = TRUE} will display
    Procrustes overlay plots from each comparison. This is the only
    step performed if input data (\code{comm}) were dissimilarities. 

    \item Scaling of the results: \code{metaMDS} will run
    \code{postMDS} for the final result. Function \code{postMDS}
    provides the following ways of ``fixing'' the indeterminacy of
    scaling and orientation of axes in NMDS: Centring moves the origin
    to the average of the axes.  Principal components rotate the
    configuration so that the variance of points is maximized on first
    dimension (with function \code{metaMDSrotate} you can
    alternatively rotate the configuration so that the first axis is
    parallel to an environmental variable).  Half-change scaling
    scales the configuration so that one unit means halving of
    community similarity from replicate similarity.  Half-change
    scaling is based on closer dissimilarities where the relation
    between ordination distance and community dissimilarity is rather
    linear; the limit is controlled by parameter \code{threshold}.  If
    there are enough points below this threshold (controlled by the
    parameter \code{nthreshold}), dissimilarities are regressed on
    distances.  The intercept of this regression is taken as the
    replicate dissimilarity, and half-change is the distance where
    similarity halves according to linear regression.  Obviously the
    method is applicable only for dissimilarity indices scaled to
    \eqn{0 \ldots 1}, such as Kulczynski, Bray-Curtis and Canberra
    indices. If half-change scaling is not used, the ordination is
    scaled to the same range as the original dissimilarities.

    \item Species scores: Function adds the species scores to the final
    solution as weighted averages using function
    \code{\link{wascores}} with given value of parameter \code{expand}.
    The expansion of  weighted averages can be undone with \code{shrink
      = TRUE} in \code{plot} or \code{scores} functions, and the
    calculation of species scores can be suppressed with \code{wascores
    = FALSE}. 
  }
}

\value{
  Function \code{metaMDS} returns an object of class
  \code{metaMDS}. The final site ordination is stored in the item
  \code{points}, and species ordination in the item \code{species}. The
  other items store the information on the steps taken by the
  function. The object has \code{print}, \code{plot}, \code{points} and
  \code{text} methods.
  Functions \code{metaMDSdist} and \code{metaMDSredist} return
  \code{\link{vegdist}} objects.
  Function \code{initMDS} returns a random
  configuration which is intended to be used within
  \code{\link[MASS]{isoMDS}} only.  Functions \code{metaMDSiter} and
  \code{postMDS} returns the result of \code{\link[MASS]{isoMDS}} with
  updated configuration.
}


\references{
  Faith, D. P, Minchin, P. R. and Belbin, L. (1987).
  Compositional dissimilarity as a robust measure of ecological
  distance. \emph{Vegetatio} 69, 57--68.

  Minchin, P.R. (1987)  An evaluation of relative robustness
  of techniques for ecological ordinations. \emph{Vegetatio} 69,
  89--107.
}

\author{ Jari Oksanen }

\note{
  Function \code{metaMDS} is a simple wrapper for
  \code{\link[MASS]{isoMDS}} and some support functions.
  You can call these support functions separately for better control
  of results.  Data transformation, dissimilarities and possible
  \code{\link{stepacross}} are made in function \code{metaMDSdist}
  which returns a dissimilarity result. Iterative search (with starting
  values from \code{initMDS}) is made in \code{metaMDSiter}.
  Processing of result configuration is done in \code{postMDS}, and
  species scores added by \code{\link{wascores}}.
  If you want
  to be more certain of reaching a global solution, you can compare
  results from several independent runs. You can also continue analysis
  from previous results or from your own configuration.  Function does
  not save the used dissimilarity matrix, but \code{metaMDSredist}
  tries to reconstruct the used dissimilarities with original data
  transformation and possible \code{\link{stepacross}}.
}
\section{Warning}{
  The calculation of \code{\link{wascores}} for species was changed in
  \pkg{vegan} version 1.12-6. They are now based on the community data
  transformed similarly as in the ordination. Previously the species
  scores always were based on the original data. You can re-establish
  the old behaviour with argument \code{old.wa = TRUE}.
}
\seealso{\code{\link[MASS]{isoMDS}}, \code{\link{decostand}},
  \code{\link{wisconsin}}, 
  \code{\link{vegdist}}, \code{\link{rankindex}}, \code{\link{stepacross}}, 
  \code{\link{procrustes}}, \code{\link{wascores}},
  \code{\link{ordiplot}}.
}

\examples{
## The recommended way of running NMDS (Minchin 1987)
##
data(dune)
library(MASS) ## isoMDS
# NMDS
sol <- metaMDS(dune)
sol
plot(sol, type="t")

## Start from previous best solution
sol2 <- metaMDS(dune, previous.best = sol)
}
\keyword{ multivariate }

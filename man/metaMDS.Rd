\encoding{UTF-8}
\name{metaMDS}
\alias{metaMDS}
\alias{metaMDSdist}
\alias{metaMDSiter}
\alias{metaMDSredist}
\alias{initMDS}
\alias{postMDS}
\alias{plot.metaMDS}
\alias{points.metaMDS}
\alias{text.metaMDS}
\alias{scores.metaMDS}

\title{Nonmetric Multidimensional Scaling with Stable Solution from
  Random Starts, Axis Scaling and Species Scores}

\description{Function \code{metaMDS} performs Nonmetric
  Multidimensional Scaling (NMDS), and tries to find a stable solution
  using several random starts. In addition, it standardizes the
  scaling in the result, so that the configurations are easier to
  interpret, and adds species scores to the site ordination. The
  \code{metaMDS} function does not provide actual NMDS, but it calls
  another function for the purpose. Currently \code{\link{monoMDS}} is
  the default choice, and it is also possible to call the
  \code{\link[MASS]{isoMDS}} (\pkg{MASS} package). }

\usage{
metaMDS(comm, distance = "bray", k = 2, try = 20, trymax = 20, 
    engine = c("monoMDS", "isoMDS"), autotransform =TRUE,
    noshare = (engine == "isoMDS"), wascores = TRUE, expand = TRUE, 
    trace = 1, plot = FALSE, previous.best,  ...)
\method{plot}{metaMDS}(x, display = c("sites", "species"), choices = c(1, 2),
    type = "p", shrink = FALSE,  ...)
\method{points}{metaMDS}(x, display = c("sites", "species"),
    choices = c(1,2), shrink = FALSE, select, ...)
\method{text}{metaMDS}(x, display = c("sites", "species"), labels, 
    choices = c(1,2), shrink = FALSE, select, ...)
\method{scores}{metaMDS}(x, display = c("sites", "species"), shrink = FALSE, 
    choices, ...)
metaMDSdist(comm, distance = "bray", autotransform = TRUE, 
    noshare = TRUE, trace = 1, commname, zerodist = "ignore", 
    distfun = vegdist, ...)
metaMDSiter(dist, k = 2, try = 20, trymax = 20, trace = 1, plot = FALSE, 
    previous.best, engine = "monoMDS", maxit = 200,
    parallel = getOption("mc.cores"), ...)   
initMDS(x, k=2)
postMDS(X, dist, pc=TRUE, center=TRUE, halfchange, threshold=0.8,
    nthreshold=10, plot=FALSE, ...)
metaMDSredist(object, ...)
}

\arguments{
  \item{comm}{Community data. Alternatively, dissimilarities either as
    a \code{\link{dist}} structure or as a symmetric square matrix. 
    In the latter case all other stages are skipped except random 
    starts and centring and pc rotation of axes. }
  \item{distance}{Dissimilarity index used in \code{\link{vegdist}}.}
  \item{k}{Number of dimensions.  NB., the number of points \eqn{n}
    should be \eqn{n > 2k + 1}{n > 2*k + 1}, and preferably higher in
    non-metric MDS.}
  \item{try, trymax}{Minimum and maximum numbers of random starts in
    search of stable solution. After \code{try} has been reached, the
    iteration will stop when two convergent solutions were found or
    \code{trymax} was reached.}
  \item{engine}{The function used for MDS. The default is to use the
    \code{\link{monoMDS}} function in \pkg{vegan}, but for backward
    compatibility it is also possible to use \code{\link{isoMDS}} of
    \pkg{MASS}.}

  \item{autotransform}{Use simple heuristics for possible data
    transformation of typical community data (see below). If you do
    not have community data, you should probably set
    \code{autotransform = FALSE}.}

  \item{noshare}{Triggering of calculation step-across or extended
    dissimilarities with function \code{\link{stepacross}}. The
    argument can be logical or a numerical value greater than zero
    and less than one. If \code{TRUE}, extended dissimilarities are
    used always when there are no shared species between some sites,
    if \code{FALSE}, they are never used. If \code{noshare} is a
    numerical value, \code{\link{stepacross}} is used when the
    proportion of site pairs with no shared species exceeds
    \code{noshare}. The number of pairs with no shared species is
    found with \code{\link{no.shared}} function, and \code{noshare}
    has no effect if input data were dissimilarities instead of
    community data.}

  \item{wascores}{Calculate species scores using function
    \code{\link{wascores}}.}
  \item{expand}{Expand weighted averages of species in
    \code{\link{wascores}}.}
  \item{trace}{Trace the function; \code{trace = 2} or higher will be
    more voluminous.}
  \item{plot}{Graphical tracing: plot interim results. You may want to set
    \code{par(ask = TRUE)} with this option.}
  \item{previous.best}{Start searches from a previous solution.}

  \item{x}{\code{metaMDS} result (or a dissimilarity structure for
  \code{initMDS}).}

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

 \item{zerodist}{Handling of zero dissimilarities: either
   \code{"fail"} or \code{"add"} a small positive value, or
   \code{"ignore"}. \code{\link{monoMDS}} accepts zero dissimilarities
   and the default is \code{zerodist = "ignore"}, but with
   \code{\link[MASS]{isoMDS}} you may need to set \code{zerodist = "add"}.}

 \item{distfun}{Dissimilarity function. Any function returning a
   \code{dist} object and accepting argument \code{method} can be used
   (but some extra arguments may cause name conflicts).}

 \item{maxit}{Maximum number of iterations in the single NMDS run;
   passed to the \code{engine} function \code{\link{monoMDS}} or
   \code{\link[MASS]{isoMDS}}.}

 \item{parallel}{Number of parallel processes or a predefined socket
   cluster.  If you use pre-defined socket clusters (say,
   \code{clus}), you must issue \code{clusterEvalQ(clus,
   library(vegan))} to make available internal \pkg{vegan}
   functions. With \code{parallel = 1} uses ordinary, non-parallel
   processing. The parallel processing is done with \pkg{parallel}
   package.}

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

 \item{\dots}{Other parameters passed to functions. Function
   \code{metaMDS} passes all arguments to its component functions
   \code{metaMDSdist}, \code{metaMDSiter}, \code{postMDS}, and to
   \code{distfun} and \code{engine}.}  }

\details{ Non-metric Multidimensional Scaling (NMDS) is commonly
  regarded as the most robust unconstrained ordination method in
  community ecology (Minchin 1987).  Function \code{metaMDS} is a
  wrapper function that calls several other functions to combine
  Minchin's (1987) recommendations into one command. The complete
  steps in \code{metaMDS} are:
  
\enumerate{ 

  \item Transformation: If the data values are larger than common
    abundance class scales, the function performs a Wisconsin double
    standardization (\code{\link{wisconsin}}).  If the values look
    very large, the function also performs \code{\link{sqrt}}
    transformation. Both of these standardizations are generally found
    to improve the results. However, the limits are completely
    arbitrary (at present, data maximum 50 triggers \code{\link{sqrt}}
    and \eqn{>9} triggers \code{\link{wisconsin}}). If you want to
    have a full control of the analysis, you should set
    \code{autotransform = FALSE} and standardize and transform data
    independently. The \code{autotransform} is intended for community
    data, and for other data types, you should set \code{autotransform
    = FALSE}. This step is perfomed using \code{metaMDSdist}.

  \item Choice of dissimilarity: For a good result, you should use
    dissimilarity indices that have a good rank order relation to
    ordering sites along gradients (Faith et al. 1987).  The default
    is Bray-Curtis dissimilarity, because it often is the test
    winner. However, any other dissimilarity index in
    \code{\link{vegdist}} can be used. Function
    \code{\link{rankindex}} can be used for finding the test winner
    for you data and gradients. The default choice may be bad if you
    analyse other than community data, and you should probably select
    an appropriate index using argument \code{distance}.  This step is
    performed using \code{metaMDSdist}.

  \item Step-across dissimilarities: Ordination may be very difficult
    if a large proportion of sites have no shared species. In this
    case, the results may be improved with \code{\link{stepacross}}
    dissimilarities, or flexible shortest paths among all sites.  The
    default NMDS \code{engine} is \code{\link{monoMDS}} which is able
    to break tied values at the maximum dissimilarity, and this often
    is sufficient to handle cases with no shared species, and
    therefore the default is not to use \code{\link{stepacross}} with
    \code{\link{monoMDS}}.  Function \code{\link[MASS]{isoMDS}} does
    not handle tied values adequately, and therefore the default is to
    use \code{\link{stepacross}} always when there are sites with no
    shared species with \code{engine = "isoMDS"}. The
    \code{\link{stepacross}} is triggered by option \code{noshare}. If
    you do not like manipulation of original distances, you should set
    \code{noshare = FALSE}.  This step is skipped if input data were
    dissimilarities instead of community data. This step is performed
    using \code{metaMDSdist}.

  \item NMDS with random starts: NMDS easily gets trapped into local
    optima, and you must start NMDS several times from random starts
    to be confident that you have found the global solution. The
    strategy in \code{metaMDS} is to first run NMDS starting with the
    metric scaling (\code{\link{cmdscale}} which usually finds a good
    solution but often close to a local optimum), or use the
    \code{previous.best} solution if supplied, and take its solution
    as the standard (\code{Run 0}). Then \code{metaMDS} starts NMDS
    from several random starts (minimum number is given by \code{try}
    and maximum number by \code{trymax}). These random starts are
    generated by \code{initMDS}. If a solution is better (has a lower
    stress) than the previous standard, it is taken as the new
    standard. If the solution is better or close to a standard,
    \code{metaMDS} compares two solutions using Procrustes analysis
    (function \code{\link{procrustes}} with option
    \code{symmetric = TRUE}). If the solutions are very similar in their
    Procrustes \code{rmse} and the largest residual is very small, the
    solutions are regarded as convergent and the better one is taken as
    the new standard.  The conditions are stringent, and you may have
    found good and relatively stable solutions although the function is
    not yet satisfied. Setting \code{trace = TRUE} will monitor the final
    stresses, and \code{plot = TRUE} will display Procrustes overlay
    plots from each comparison. This step is performed using
    \code{metaMDSiter}. This is the only step performed if input data
    (\code{comm}) were dissimilarities.

  \item Scaling of the results: \code{metaMDS} will run \code{postMDS}
    for the final result. Function \code{postMDS} provides the
    following ways of \dQuote{fixing} the indeterminacy of scaling and
    orientation of axes in NMDS: Centring moves the origin to the
    average of the axes; Principal components rotate the configuration
    so that the variance of points is maximized on first dimension
    (with function \code{\link{MDSrotate}} you can alternatively rotate
    the configuration so that the first axis is parallel to an
    environmental variable); Half-change scaling scales the
    configuration so that one unit means halving of community
    similarity from replicate similarity.  Half-change scaling is
    based on closer dissimilarities where the relation between
    ordination distance and community dissimilarity is rather linear
    (the limit is set by argument \code{threshold}). If there are
    enough points below this threshold (controlled by the parameter
    \code{nthreshold}), dissimilarities are regressed on distances.
    The intercept of this regression is taken as the replicate
    dissimilarity, and half-change is the distance where similarity
    halves according to linear regression.  Obviously the method is
    applicable only for dissimilarity indices scaled to \eqn{0 \ldots
    1}, such as Kulczynski, Bray-Curtis and Canberra indices. If
    half-change scaling is not used, the ordination is scaled to the
    same range as the original dissimilarities.

  \item Species scores: Function adds the species scores to the final
    solution as weighted averages using function
    \code{\link{wascores}} with given value of parameter
    \code{expand}.  The expansion of weighted averages can be undone
    with \code{shrink = TRUE} in \code{plot} or \code{scores}
    functions, and the calculation of species scores can be suppressed
    with \code{wascores = FALSE}.
} 
}

\section{Convergence Problems}{

    The function tries hard to find two convergent solutions, but it
    may fail. With default \code{engine = "monoMDS"} the function will
    tabulate the stopping criteria used, so that you can see which
    criterion should be made more stringent. The criteria can be given
    as arguments to \code{metaMDS} and their current values are
    described in \code{\link{monoMDS}}. In particular, if you reach
    the maximum number of iterations, you should increase the value of
    \code{maxit}. You may ask for a larger number of random starts
    without losing the old ones giving the previous solution in
    argument \code{previous.best}.

    In addition to too slack convergence criteria and too low number
    of random starts, wrong number of dimensions (argument \code{k})
    is the most common reason for not finding convergent
    solutions. NMDS is usually run with a low number dimensions
    (\code{k=2} or \code{k=3}), and for complex data increasing
    \code{k} by one may help. If you run NMDS with much higher number
    of dimensions (say, \code{k=10} or more), you should reconsider
    what you are doing and drastically reduce \code{k}. For very
    heterogeneous data sets with partial disjunctions, it may help to
    set \code{stepacross}, but for most data sets the default
    \code{weakties = TRUE} is sufficient.

    Please note that you can give all arguments of other
    \code{metaMDS*} functions and NMDS engine (default
    \code{\link{monoMDS}}) in your \code{metaMDS} command,and you
    should check documentation of these functions for details.

}

\value{ Function \code{metaMDS} returns an object of class
  \code{metaMDS}. The final site ordination is stored in the item
  \code{points}, and species ordination in the item \code{species},
  and the stress in item \code{stress} (NB, the scaling of the stress
  depends on the \code{engine}: \code{\link[MASS]{isoMDS}} uses
  percents, and \code{\link{monoMDS}} proportions in the range \eqn{0
  \ldots 1}). The other items store the information on the steps taken
  and the items returned by the \code{engine} function. The object has
  \code{print}, \code{plot}, \code{points} and \code{text} methods.
  Functions \code{metaMDSdist} and \code{metaMDSredist} return
  \code{\link{vegdist}} objects.  Function \code{initMDS} returns a
  random configuration which is intended to be used within
  \code{\link[MASS]{isoMDS}} only.  Functions \code{metaMDSiter} and
  \code{postMDS} returns the result of NMDS with updated
  configuration.  }


\references{
  Faith, D. P, Minchin, P. R. and Belbin, L. (1987).
  Compositional dissimilarity as a robust measure of ecological
  distance. \emph{Vegetatio} 69, 57--68.

  Minchin, P.R. (1987)  An evaluation of relative robustness
  of techniques for ecological ordinations. \emph{Vegetatio} 69,
  89--107.
}

\author{ Jari Oksanen }

\note{ Function \code{metaMDS} is a simple wrapper for an NMDS engine
  (either \code{\link{monoMDS}} or \code{\link[MASS]{isoMDS}}) and
  some support functions (\code{metaMDSdist},
  \code{\link{stepacross}}, \code{metaMDSiter}, \code{initMDS},
  \code{postMDS}, \code{\link{wascores}}).  You can call these support
  functions separately for better control of results.  Data
  transformation, dissimilarities and possible
  \code{\link{stepacross}} are made in function \code{metaMDSdist}
  which returns a dissimilarity result. Iterative search (with
  starting values from \code{initMDS} with \code{\link{monoMDS}}) is
  made in \code{metaMDSiter}.  Processing of result configuration is
  done in \code{postMDS}, and species scores added by
  \code{\link{wascores}}.  If you want to be more certain of reaching
  a global solution, you can compare results from several independent
  runs. You can also continue analysis from previous results or from
  your own configuration.  Function may not save the used
  dissimilarity matrix (\code{\link{monoMDS}} does), but
  \code{metaMDSredist} tries to reconstruct the used dissimilarities
  with original data transformation and possible
  \code{\link{stepacross}}.

  The \code{metaMDS} function was designed to be used with community
  data.  If you have other type of data, you should probably set some
  arguments to non-default values: probably at least \code{wascores},
  \code{autotransform} and \code{noshare} should be \code{FALSE}. If
  you have negative data entries, \code{metaMDS} will set the previous
  to \code{FALSE} with a warning.  
}

\section{Warning}{\code{metaMDS} uses \code{\link{monoMDS}} as its
  NMDS \code{engine} from \pkg{vegan} version 2.0-0, when it replaced
  the \code{\link[MASS]{isoMDS}} function. You can set argument
  \code{engine} to select the old engine.}

\seealso{\code{\link{monoMDS}} (and \code{\link[MASS]{isoMDS}}), 
  \code{\link{decostand}},
  \code{\link{wisconsin}}, 
  \code{\link{vegdist}}, \code{\link{rankindex}}, \code{\link{stepacross}}, 
  \code{\link{procrustes}}, \code{\link{wascores}}, \code{\link{MDSrotate}},
  \code{\link{ordiplot}}.
}

\examples{
## The recommended way of running NMDS (Minchin 1987)
##
data(dune)
# Global NMDS using monoMDS
sol <- metaMDS(dune)
sol
plot(sol, type="t")
## Start from previous best solution
sol <- metaMDS(dune, previous.best = sol)
## Local NMDS and stress 2 of monoMDS
sol2 <- metaMDS(dune, model = "local", stress=2)
sol2
## Use Arrhenius exponent 'z' as a binary dissimilarity measure
sol <- metaMDS(dune, distfun = betadiver, distance = "z")
sol
}
\keyword{ multivariate }

\name{mrpp}
\alias{mrpp}
\alias{print.mrpp}

\title{ Multi Response Permutation Procedure of Within- versus
  Among-Group Dissimilarities}

\description{ Multiple Response Permutation Procedure (MRPP) provides a
test of whether there is a significant difference between two or more
groups of sampling units.  }

\usage{
mrpp(dat, grouping, permutations = 1000, distance = "euclidean",
     weight.type = 1, strata)
}

\arguments{
  \item{dat}{data matrix or data frame in which rows are samples and
    columns are response variable(s), or a dissimilarity object or a
    symmetric square matrix of dissimilarities.} 
  \item{grouping}{ Factor or numeric index for grouping observations.}
  \item{permutations}{Number of permutations to assess the significance
    of the MRPP statistic, \eqn{delta}.} 
  \item{distance}{Choice of distance metric that measures the
    dissimilarity between two observations . See \code{\link{vegdist}} for
    options.  This will be used if \code{dat} was not a dissimilarity
    structure of a symmetric square matrix.}  
  \item{weight.type}{ choice of group weights. See Details below for options.}
  \item{strata}{An integer vector or factor specifying the strata for
    permutation. If supplied, observations are permuted only within the
    specified strata.}
}

\details{ Multiple Response Permutation Procedure (MRPP) provides a test
of whether there is a significant difference between two or more groups
of sampling units. This difference may be one of location (differences
in mean) or one of spread (differences in within-group
distance). Function \code{mrpp} operates on a \code{data.frame} matrix
where rows are observations and responses data matrix. The response(s)
may be uni- or multivariate. The method is philosophically and
mathematically allied with analysis of variance, in that it compares
dissimilarities within and among groups. If two groups of sampling units
are really different (e.g. in their species composition), then average
of the within-group compositional dissimilarities ought to be less than
the average of the dissimilarities between two random collection of
sampling units drawn from the entire population. 

The mrpp statistic \eqn{\delta} is simply the overall weighted mean of
within-group means of the pairwise dissimilarities among sampling
units. The correct choice of group weights is currently not clear. The
\code{mrpp} function offers three choices: (1) group size (\eqn{n}), (2) a
degrees-of-freedom analogue (\eqn{n-1}), and (3) a weight that is the number
of unique distances calculated among \eqn{n} sampling units (\eqn{n(n-1)/2}).

The \code{mrpp} algorithm first calculates all pairwise distances in the
entire dataset, then calculates \eqn{\delta}. It then permutes the
sampling units and their associated pairwise distances, and recalculates
a \eqn{\delta} based on the permuted data. It repeats the permutation
step \code{permutations} times. The significance test is simply the
fraction of permuted deltas that are less than the observed delta, with
a small sample correction. The function also calculates the
change-corrected within-group agreement
\eqn{A = 1 -\delta/E(\delta)}, where \eqn{E(\delta)} is the expected
\eqn{\delta} assessed as the average of permutations.

With \code{weight.type = 3}, the function also calculates classification
strength (Van Sickle 1997) which is defined as the difference between
average between group dissimilarities and within group
dissimilarities. With \code{weight.type = 3} the classification strength
is a simple transformation of \eqn{\delta}, and has the same permutation
significance.

If the first argument \code{dat} can be interpreted as dissimilarities,
they will be used directly. In other cases the function treats
\code{dat} as observations, and uses \code{\link{vegdist}} to find 
the dissimilarities.  The default \code{distance} is Euclidean as in the
traditional use of the method, but other dissimilarities in
\code{\link{vegdist}} also are available. }

\value{
The function returns a list of class mrpp with following items:
  \item{call }{	Function call.}
  \item{delta }{The overall weighted mean of group mean distances.}
  \item{E.delta}{expected delta, under the null hypothesis of no group
    structure. This is the mean of original dissimilarities.}
  \item{CS}{Classification strength (Van Sickle 1997) with
    \code{weight.type = 3} and \code{NA} with other weights.}
  \item{n}{Number of observations in each class.}
  \item{classdelta}{Mean dissimilarities within classes. The overall
    \eqn{\delta} is the weighted average of these values with given
    \code{weight.type}}. 
  \item{Pvalue}{Significance of the test.}
  \item{A}{A chance-corrected estimate of the proportion of the distances
    explained by group identity; a value analogous to a coefficient of 
    determination in a linear model.} 
  \item{distance}{Choice of distance metric used; the "method" entry of
    the dist object.}
  \item{weight.type}{The choice of group weights used.}
  \item{boot.deltas}{The vector of "permuted deltas," the deltas
    calculated from each of the permuted datasets.}
  \item{permutations}{The number of permutations used.}
}
\references{
  B. McCune and J. B. Grace. 2002. \emph{Analysis of Ecological
  Communities.} MjM  Software Design, Gleneden Beach, Oregon, USA.


  P. W. Mielke and K. J. Berry. 2001. \emph{Permutation Methods: A
  Distance  Function Approach.} Springer Series in
  Statistics. Springer.  

  J. Van Sickle 1997. Using mean similarity dendrograms to evaluate
  classifications. \emph{Journal of Agricultural, Biological, and Environmental
  Statistics} 2:370-388.

}
\author{
M. Herny H. Stevens \email{HStevens@muohio.edu} and Jari Oksanen.
}
\note{ 
This difference may be one of location (differences in mean) or one of
spread (differences in within-group distance). That is, it may find a
significant difference between two groups simply because one of those
groups has a greater dissimilarities among its sampling units. Most
\code{mrpp} models can be analysed with \code{\link{adonis}} which seems
not suffer from the same problems as \code{mrpp} and is a more robust
alternative. 
}
\seealso{
  \code{\link{anosim}} for a similar test based on ranks, and
  \code{\link{mantel}} for comparing dissimilarities against continuous
  variables, and
  \code{\link{vegdist}} for obtaining dissimilarities,
  \code{\link{adonis}} is a more robust alternative in most cases.
}
\examples{
data(dune)
data(dune.env)
dune.mrpp <- mrpp(dune, dune.env$Management)
dune.mrpp

# Save and change plotting parameters
def.par <- par(no.readonly = TRUE)
layout(matrix(1:2,nr=1))

plot(dune.ord <- metaMDS(dune), type="text", display="sites" )
ordihull(dune.ord, dune.env$Management)

with(dune.mrpp, {
  fig.dist <- hist(boot.deltas, xlim=range(c(delta,boot.deltas)), 
                 main="Test of Differences Among Groups")
  abline(v=delta); 
  text(delta, 2*mean(fig.dist$counts), adj = -0.5,
     expression(bold(delta)), cex=1.5 )  }
)
par(def.par)
}
\keyword{ multivariate }
\keyword{ nonparametric }
\keyword{htest}

\name{sppscores}
\alias{sppscores}
\alias{sppscores<-}
\alias{sppscores<-.dbrda}
\alias{sppscores<-.capscale}
\alias{sppscores<-.metaMDS}
\alias{sppscores<-.monoMDS}
\alias{sppscores<-.wcmdscale}

\title{
Add or Replace Species Scores in Distance-Based Ordination
}

\description{

  Distance-based ordination (\code{\link{dbrda}},
  \code{\link{capscale}}, \code{\link{metaMDS}}, \code{\link{monoMDS}},
  \code{\link{wcmdscale}}) has no information on species, but some
  methods may add species scores if community data were
  available. However, the species scores may be missing (and they always
  are in \code{\link{dbrda}} and \code{\link{wcmdscale}}), or they may
  not have a close relation to used dissimilarity index. This function
  will add the species scores or replace the existing species scores in
  distance-based methods.

}

\usage{
sppscores(object) <- value
}

\arguments{
  \item{object}{Ordination result from \code{\link{capscale}},
    \code{\link{dbrda}}, \code{\link{metaMDS}}, \code{\link{monoMDS}},
    \code{\link{pco}} or \code{\link{wcmdscale}}.}
  \item{value}{Community data to find the species scores.}
}

\details{

  Distances have no information on species (columns, variables), and
  hence distance-based ordination has no information on species
  scores. However, the species scores can be added as supplementary
  information after the analysis to help the interpretation of
  results. Some ordination methods (\code{\link{capscale}},
  \code{\link{metaMDS}}) can supplement the species scores during the
  analysis if community data were available in the analysis.

  In \code{\link{capscale}} the species scores are found by projecting
  the community data to site ordination (linear combination scores),
  and the scores are accurate if the analysis used Euclidean
  distances. If the dissimilarity index can be expressed as Euclidean
  distances of transformed data (for instance, Chord and Hellinger
  Distances), the species scores based on transformed data will be
  accurate, but the function still finds the dissimilarities with
  untransformed data. Usually community dissimilarities differ in two
  significant ways from Euclidean distances: They are bound to maximum
  1, and they use absolute differences instead of squared
  differences. In such cases, it may be better to use species scores
  that are transformed so that their Euclidean distances have a good
  linear relation to used dissimilarities. It is often useful to
  standardize data so that each row has unit total, and perform
  squareroot transformation to damp down the effect of squared
  differences (see Examples).

  Functions \code{\link{dbrda}} and \code{\link{wcmdscale}} never find
  the species scores, but they mathematically similar to
  \code{\link{capscale}}, and similar rules should be followed when
  supplementing the species scores.

  Functions for species scores in \code{\link{metaMDS}} and
  \code{\link{monoMDS}} use weighted averages (\code{\link{wascores}})
  to find the species scores. These have better relationship with most
  dissimilarities than the projection scores used in metric ordination,
  but similar transformation of the community data should be used both
  in dissimilarities and in species scores.

}

\value{

  Replacement function adds the species scores or replaces the old
  scores in the ordination object.

}

\author{
  Jari Oksanen
}

\seealso{

  Function \code{\link{envfit}} finds similar scores, but based on
  correlations. The species scores for non-metric ordination use
  \code{\link{wascores}} which can also used directly on any
  ordination result.
 
}

\examples{
data(BCI, BCI.env)
mod <- dbrda(vegdist(BCI) ~ Habitat, BCI.env)
## add species scores
sppscores(mod) <- BCI
## Euclidean distances of BCI differ from used dissimilarity
plot(vegdist(BCI), dist(BCI))
## more linear relationship
plot(vegdist(BCI), dist(sqrt(decostand(BCI, "total"))))
## better species scores
sppscores(mod) <- sqrt(decostand(BCI, "total"))
}

\keyword{ multivariate }


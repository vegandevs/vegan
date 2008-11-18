\name{vegemite}
\alias{vegemite}
\alias{coverscale}

\title{Prints a Compact, Ordered Vegetation Table }
\description{
  The function prints a compact vegetation table, where species are
  rows, and each site takes only one column without spaces.  The
  vegetation table can be ordered by explicit indexing, by environmental
  variables or results from an ordination or cluster analysis.
}
\usage{
vegemite(x, use, scale, sp.ind, site.ind, zero=".", select ,...)
coverscale(x, scale=c("Braun.Blanquet", "Domin", "Hult", "Hill", "fix","log"),
           maxabund)
}

\arguments{
  \item{x}{Vegetation data. }
  \item{use}{Either a vector, or an object from \code{cca},
    \code{decorana} \emph{etc.} or \code{hclust} or a
    \code{\link{dendrogram}} for ordering sites and species.} 
  \item{sp.ind}{Species indices. }
  \item{site.ind}{Site indices. }
  \item{zero}{Character used for zeros. }
  \item{select}{Select a subset of sites.  This can be a logical vector
    (\code{TRUE} for selected sites), or a vector of indices of selected
    sites.  The order of indices does not influence results, but you
    must specify \code{use} or \code{site.ind} to reorder sites.
    }
  \item{scale}{Cover scale used (can be abbreviated).}
  \item{maxabund}{Maximum abundance used with \code{scale = "log"}.
    Data maximum in the \code{select}ed subset will be usd if this is
    missing.}
  \item{...}{Arguments passed to \code{coverscale} (i.e., \code{maxabund}).}
}
\details{
  The function prints a traditional vegetation table.
  Unlike in ordinary data matrices, species are used as rows and sites
  as columns.  The table is printed in compact form:  only one character
  can be used for abundance, and there are no spaces between
  columns. Species with no occurrences are dropped from the table.

  The parameter \code{use} can be a vector or an object from
  \code{\link{hclust}}, a \code{\link{dendrogram}} or any ordination
  result recognized by \code{\link{scores}} (all ordination methods in
  \pkg{vegan} and some of those not in \pkg{vegan}).
  If \code{use} is a vector, it is used
  for ordering sites.  If \code{use} is an object from ordination, both
  sites and species are arranged by the first axis.  
  When \code{use} is an
  object from \code{\link{hclust}} or a \code{\link{dendrogram}}, the
  sites are ordered similarly
  as in the cluster dendrogram.
  If ordination methods provide species scores, these are used for
  ordering species.  In all cases where species scores are missing,
  species are ordered by their weighted averages (\code{\link{wascores}})
  on site scores. There is no natural, unique ordering in hierarchic
  clustering, but in some cases species are still nicely ordered (please
  note that you can \code{\link{reorder.dendrogram}} to have such a
  natural order).
  Alternatively, species and sites can be ordered explicitly giving
  their indices or names in parameters \code{sp.ind} and
  \code{site.ind}.  If these are given, they take precedence over
  \code{use}. A subset of sites can be displayed using argument
  \code{select}, but this cannot be used to order sites, but you still
  must give \code{use} or \code{site.ind}.

  If \code{scale} is given, \code{vegemite} calls
  \code{coverscale} to transform percent cover
  scale or some other scales into traditional class scales used in
  vegetation science (\code{coverscale} can be called directly, too).
  Braun-Blanquet and Domin scales are actually not
  strict cover scales, and the limits used for codes \code{r} and
  \code{+} are arbitrary.  Scale \code{Hill} may be
  inappropriately named, since Mark O. Hill probably never intended this
  as a cover scale.  However, it is used as default 'cut levels' in his
  \code{TWINSPAN}, and surprisingly many users stick to this default,
  and this is a \emph{de facto} standard in publications.  All
  traditional
  scales assume that values are cover percentages with maximum 100.
  However, non-traditional alternative \code{log} can be used with any
  scale range.  Its class limits are integer powers of 1/2 of the
  maximum (argument \code{maxabund}), with \code{+} used for non-zero entries
  less than 1/512 of the maximum (\code{log} stands alternatively for
  logarithmic or logical).  Scale \code{fix} is intended for `fixing'
  10-point scales: it truncates scale values to integers, and replaces
  10 with \code{X} and positive values below 1 with \code{+}. 
}
\value{
  The function is used mainly to print a table, but it returns
  (invisibly) a list
  with items:
  \item{species}{Ordered species indices}
  \item{sites}{Ordered site indices}
  These items can be used as arguments \code{sp.ind} and \code{site.ind}
  to reproduce the table. In addition to the proper table, the function
  prints the numbers of species and sites and the name of the used cover
  scale at the end.
}
\references{ The cover scales are presented in many textbooks of vegetation
  science; I used:

  Shimwell, D.W. (1971) \emph{The Description and Classification of
  Vegetation}. Sidgwick & Jackson.
}
\author{Jari Oksanen}

\seealso{\code{\link{cut}} and \code{\link{approx}} for making your own
  `cover scales', \code{\link{wascores}} for weighted averages.
  }

\note{ This function was called \code{vegetab} in older versions of
  \code{vegan}.  The new name was chosen  because the output is so
  compact (and to avoid confusion with the \code{vegtab} function in the
  \pkg{labdsv} package).
    }
\examples{
data(varespec)
## Print only more common species 
freq <- apply(varespec > 0, 2, sum)
vegemite(varespec, scale="Hult", sp.ind = freq > 10)
## Order by correspondence analysis, use Hill scaling and layout:
dca <- decorana(varespec)
vegemite(varespec, dca, "Hill", zero="-")
## Show one class from cluster analysis, but retain the ordering above
clus <- hclust(vegdist(varespec))
cl <- cutree(clus, 3)
sel <- vegemite(varespec, use=dca, select = cl == 3, scale="Br")
# Re-create previous
vegemite(varespec, sp=sel$sp, site=sel$site, scale="Hult")
}
\keyword{ print }
\keyword{ manip }


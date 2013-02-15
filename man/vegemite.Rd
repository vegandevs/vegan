\name{vegemite}
\alias{vegemite}
\alias{tabasco}
\alias{coverscale}

\title{Display Compact Ordered Community Tables }
\description{
  The functions display compact community tables.  Function
  \code{vegemite} prints text tables where species are rows, and each
  site takes only one column without spaces.  Function \code{tabasco}
  provides interface for \code{\link{heatmap}} for a colour map of
  the data. The community table can be ordered by explicit indexing,
  by environmental variables or results from an ordination or cluster
  analysis.
}
\usage{
vegemite(x, use, scale, sp.ind, site.ind, zero=".", select ,...)
tabasco(x, use, sp.ind = NULL, site.ind = NULL, select,
    Rowv = TRUE, Colv = TRUE, ...) 
coverscale(x, scale=c("Braun.Blanquet", "Domin", "Hult", "Hill", "fix","log"),
           maxabund)
}

\arguments{
  \item{x}{Vegetation data. }
  \item{use}{Either a vector, or an object from \code{cca},
    \code{decorana} \emph{etc.} or \code{hclust} or a
    \code{\link{dendrogram}} for ordering sites and species.} 
  \item{sp.ind}{Species indices. In \code{tabasco}, this can also be an
    \code{\link{hclust}} tree or a \code{\link{dendrogram}}. }
  \item{site.ind}{Site indices. }
  \item{zero}{Character used for zeros. }
  \item{select}{Select a subset of sites.  This can be a logical vector
    (\code{TRUE} for selected sites), or a vector of indices of selected
    sites.  The order of indices does not influence results, but you
    must specify \code{use} or \code{site.ind} to reorder sites.
  }
  \item{Rowv, Colv}{Reorder dendrograms for the rows (sites) or columns
    (species) of \code{x}.  If the arguments are \code{TRUE}, row
    dendrograms are ordered by the first axis of correspondence
    analysis, and column dendrograms by the weighted average
    (\code{\link{wascores}}) of the order in the row.  Alternatively
    these can be vectors that are used to reorder the dendrogram. }
  \item{scale}{Cover scale used (can be abbreviated).}
  \item{maxabund}{Maximum abundance used with \code{scale = "log"}.
    Data maximum in the \code{select}ed subset will be used if this is
    missing.}
  \item{...}{Arguments passed to \code{coverscale} (i.e., \code{maxabund}) in
    \code{vegemite} and to \code{\link{heatmap}} in \code{tabasco}.}
}
\details{
  The function \code{vegemite} prints a traditional community table.
  Unlike in ordinary data matrices, species are used as rows and sites
  as columns.  The table is printed in compact form:  only one character
  can be used for abundance, and there are no spaces between
  columns. Species with no occurrences are dropped from the table.

  Function \code{tabasco} produces a similar table as \code{vegemite}
  using \code{\link{heatmap}}, where abundances are coded by heatmap
  colours.

  The parameter \code{use} can be a vector or an object from
  \code{\link{hclust}}, a \code{\link{dendrogram}} or any ordination
  result recognized by \code{\link{scores}} (all ordination methods in
  \pkg{vegan} and some of those not in \pkg{vegan}). The
  \code{\link{hclust}} and \code{\link{dendrogram}} must be for
  sites. The dendrogram is displayed above the sites in
  \code{tabasco}.  No dendrogram for species is displayed, except when
  given in \code{sp.ind}.

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

  Species and sites can be ordered explicitly giving their indices or
  names in parameters \code{sp.ind} and \code{site.ind}.  If these are
  given, they take precedence over \code{use}. A subset of sites can be
  displayed using argument \code{select}, but this cannot be used to
  order sites, but you still must give \code{use} or \code{site.ind}.
  However, \code{tabasco} makes two exceptions: \code{site.ind} and
  \code{select} cannot be used when \code{use} is an
  \code{\link{hclust}} tree or a \code{\link{dendrogram}}. In addition,
  the \code{sp.ind} can also be an \code{\link{hclust}} tree or a
  \code{\link{dendrogram}}, and in that case the dendrogram is plotted
  on the left side of the \code{\link{heatmap}}.

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
  The functions are used mainly to display a table, but it returns
  (invisibly) a list
  with items:
  \item{species}{Ordered species indices}
  \item{sites}{Ordered site indices}
  These items can be used as arguments \code{sp.ind} and \code{site.ind}
  to reproduce the table. In addition to the proper table, \code{vegemite}
  prints the numbers of species and sites and the name of the used cover
  scale at the end.
}
\references{ The cover scales are presented in many textbooks of vegetation
  science; I used:

  Shimwell, D.W. (1971) \emph{The Description and Classification of
  Vegetation}. Sidgwick & Jackson.
}
\author{Jari Oksanen}

\seealso{\code{\link{cut}} and \code{\link{approx}} for making your
  own \sQuote{cover scales} for \code{vegemite}.  Function
  \code{tabasco} is based on \code{\link{heatmap}}. Both functions
  order species with weighted averages using \code{\link{wascores}}.
  }

\note{   

  The name \code{vegemite} was chosen because the output is so
  compact, and the \code{tabasco} because it is just as compact, but
  uses heat colours.
    
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

## Abundance values have such a wide range that they must be rescaled
## or all abundances will not be shown in tabasco
tabasco(decostand(varespec, "log"), clus)
## reorder dendrogram by the first DCA axis
clus <- as.dendrogram(clus)
clus <- reorder(clus, scores(dca, choices=1, display="sites"))
tabasco(decostand(varespec, "log"), clus)

}
\keyword{ print }
\keyword{ manip }


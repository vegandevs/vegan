\name{treedive}
\alias{treedive}
\alias{treeheight}
\alias{treedist}  

\title{Functional Diversity and Community Distances from Species Trees}

\description{ Functional diversity is defined as the total branch
  length in a trait dendrogram connecting all species, but excluding
  the unnecessary root segments of the tree (Petchey and Gaston
  2006). Tree distance is the increase in total branch length when
  combining two sites.  }

\usage{
treedive(comm, tree, match.force = TRUE, verbose = TRUE)
treeheight(tree)
treedist(x, tree, relative = TRUE, match.force = TRUE, ...)
}

\arguments{
  \item{comm, x}{Community data frame or matrix.}
  \item{tree}{A dendrogram which for \code{treedive} must be for species
    (columns).}
  \item{match.force}{Force matching of column names in data
    (\code{comm}, \code{x}) and labels in \code{tree}. If \code{FALSE},
    matching only happens when dimensions differ (with a warning or
    message). The order of data must match to the order in \code{tree}
    if matching by names is not done.}
  \item{verbose}{Print diagnostic messages and warnings.}
  \item{relative}{Use distances relative to the height of combined tree.}
  \item{\dots}{Other arguments passed to functions (ignored).}
}
\details{
  Function \code{treeheight} finds the sum of lengths of connecting
  segments in a dendrogram produced by \code{\link{hclust}}, or other
  dendrogram that can be coerced to a correct type using
  \code{\link{as.hclust}}. When applied to a clustering of species
  traits, this is a measure of functional diversity (Petchey and Gaston
  2002, 2006), and when applied to phylogenetic trees this is
  phylogenetic diversity.

  Function \code{treedive} finds the \code{treeheight} for each site
  (row) of a community matrix. The function uses a subset of
  dendrogram for those species that occur in each site, and excludes
  the tree root if that is not needed to connect the species (Petchey
  and Gaston 2006). The subset of the dendrogram is found by first
  calculating \code{\link{cophenetic}} distances from the input
  dendrogram, then reconstructing the dendrogram for the subset of the
  cophenetic distance matrix for species occurring in each
  site. Diversity is 0 for one species, and \code{NA} for empty
  communities.

  Function \code{treedist} finds the dissimilarities among
  trees. Pairwise dissimilarity of two trees is found by combining
  species in a common tree and seeing how much of the tree height is
  shared and how much is unique. With \code{relative = FALSE} the
  dissimilarity is defined as \eqn{2 (A \cup B) - A - B}{2AB - A - B}, where
  \eqn{A} and \eqn{B} are heights of component trees and
  \eqn{A \cup B}{AB} is the height of the combined tree. With \code{relative = TRUE}
  the dissimilarity is \eqn{(2(A \cup B)-A-B)/(A \cup B)}{(2AB - A - B)/(AB)}. 
  Although the latter formula is similar to
  Jaccard dissimilarity (see \code{\link{vegdist}},
  \code{\link{designdist}}), it is not in the range \eqn{0 \ldots 1}{0
  \dots 1}, since combined tree can add a new root. When two zero-height
  trees are combined into a tree of above zero height, the relative
  index attains its maximum value \eqn{2}. The dissimilarity is zero
  from a combined zero-height tree.

  The functions need a dendrogram of species traits or phylogenies as an
  input. If species traits contain \code{\link{factor}} or
  \code{\link{ordered}} factor variables, it is recommended to use Gower
  distances for mixed data (function \code{\link[cluster]{daisy}} in
  package \pkg{cluster}), and usually the recommended clustering method
  is UPGMA (\code{method = "average"} in function \code{\link{hclust}})
  (Podani and Schmera 2006). Phylogenetic trees can be changed into
  dendrograms using function \code{as.hclust.phylo} in the
  \pkg{ape} package.

  It is possible to analyse the non-randomness of tree diversity
  using \code{\link{oecosimu}}. This needs specifying an adequate Null
  model, and the results will change with this choice.
}

\value{ 
  A vector of diversity values or a single tree height, or a
  dissimilarity structure that inherits from \code{\link{dist}} and
  can be used similarly.  
}

\references{
  Lozupone, C. and Knight, R. 2005. UniFrac: a new phylogenetic method
  for comparing microbial communities. \emph{Applied and Environmental
  Microbiology} 71, 8228--8235.

  Petchey, O.L. and Gaston, K.J. 2002. Functional diversity (FD), species
  richness and community composition. \emph{Ecology Letters} 5,
  402--411.
  
  Petchey, O.L. and Gaston, K.J. 2006. Functional diversity: back to
  basics and looking forward. \emph{Ecology Letters} 9, 741--758.

  Podani J. and Schmera, D. 2006. On dendrogram-based methods of
  functional diversity. \emph{Oikos} 115, 179--185.
}

\author{Jari Oksanen}

\seealso{
  Function \code{treedive} is similar to the phylogenetic
  diversity function \code{pd} in the package \pkg{picante}, but
  excludes tree root if that is not needed to connect species. Function
  \code{treedist} is similar to the phylogenetic similarity
  \code{phylosor} in the package \pkg{picante}, but excludes
  unneeded tree root and returns distances instead of similarities.

  \code{\link{taxondive}} is something very similar from another bubble.
}

\examples{

## There is no data set on species properties yet, and we demonstrate
## the methods using phylogenetic trees
data(dune)
data(dune.phylodis)
cl <- hclust(dune.phylodis)
treedive(dune, cl)
## Significance test using Null model communities.
## The current choice fixes numbers of species and picks species
## proportionally to their overall frequency
oecosimu(dune, treedive, "r1", tree = cl, verbose = FALSE)
## Phylogenetically ordered community table
dtree <- treedist(dune, cl)
tabasco(dune, hclust(dtree), cl)
## Use tree distances in distance-based RDA
dbrda(dtree ~ 1)
}

\keyword{ univar }


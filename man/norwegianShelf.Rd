\name{norwegianShelf}
\alias{norwegianShelf}
\docType{data}
\title{Norwegian Continental Shelf Macrofauna Data}
\usage{
data("norwegianShelf")
}
\description{
  Counts of 809 benthic soft-sediment macrofaunal species at 101 sites
  from five areas along the Norwegian continental shelf. Anderson et al.
  (2017) used presence/absence versions of these data to illustrate the
  modified PERMANOVA statistic for heterogeneous dispersions.
}
\format{
  A list with three components:
  \describe{
    \item{comm}{An integer matrix with 101 sites as rows and 809 species
      as columns.}
    \item{area}{A factor assigning the 101 sites to five sampling areas.}
    \item{site}{Integer site identifiers from 1 to 101.}
  }
}
\details{
  The published presence/absence analysis is obtained with
  \code{1L * (norwegianShelf$comm > 0L)}.
}
\source{
  Supporting Data S1 for Anderson et al. (2017),
  \doi{10.1111/anzs.12176}.
}
\references{
  Anderson, M.J., Walsh, D.C.I., Clarke, K.R., Gorley, R.N. and Guerra-
  Castro, E. (2017). Some solutions to the multivariate Behrens--Fisher
  problem for dissimilarity-based analyses. \emph{Australian & New
  Zealand Journal of Statistics} 59, 57--79.
}
\examples{
data("norwegianShelf")
str(norwegianShelf, max.level = 1)
}
\keyword{datasets}

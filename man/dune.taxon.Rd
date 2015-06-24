\encoding{UTF-8}
\name{dune.taxon}
\alias{dune.taxon}
\alias{dune.phylodis}
\docType{data}
\title{Taxonomic Classification and Phylogeny of Dune Meadow Species}
\description{
  Classification table of the species in the \code{\link{dune}} data
  set.
}
\usage{
  data(dune.taxon)
  data(dune.phylodis)
}
\format{
  \code{dune.taxon} is data frame with 30 species (rows) classified
  into five taxonomic levels (columns). \code{dune.phylodis} is a
  \code{\link{dist}} object of estimated coalescence ages extracted
  from \url{http://datadryad.org/resource/doi:10.5061/dryad.63q27}
  (Zanne et al. 2014) using tools in packages \pkg{ape} and
  \pkg{phylobase}.
}

\details{
  The classification of vascular plants is based on APG (2009), and
  that of mosses on Hill et al. (2006).
}
\references{
  APG [Angiosperm Phylogeny Group] (2009) An update of the Angiosperm
  Phylogeny Group classification for the orders and families of flowering
  plants: APG III. \emph{Bot. J. Linnean Soc.} \strong{161}: 105--121.

  Hill, M.O et al. (2006) An annotated checklist of the mosses of Europe
  and Macaronesia. \emph{J. Bryology} \strong{28}: 198--267.

  Zanne A.E., Tank D.C., Cornwell, W.K., Eastman J.M., Smith, S.A.,
  FitzJohn, R.G., McGlinn, D.J., O’Meara, B.C., Moles, A.T., Reich,
  P.B., Royer, D.L., Soltis, D.E., Stevens, P.F., Westoby, M., Wright,
  I.J., Aarssen, L., Bertin, R.I., Calaminus, A., Govaerts, R.,
  Hemmings, F., Leishman, M.R., Oleksyn, J., Soltis, P.S., Swenson,
  N.G., Warman, L. & Beaulieu, J.M. (2014) Three keys to the radiation
  of angiosperms into freezing environments. \emph{Nature} 506, 89--92.

}

\seealso{Functions \code{\link{taxondive}}, \code{\link{treedive}},
  and \code{\link{treedist}} use these data sets. }

\examples{ 
  data(dune.taxon) 
  data(dune.phylodis)
}
\keyword{datasets}

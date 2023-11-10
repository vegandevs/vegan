\encoding{UTF-8}
\name{BCI}
\alias{BCI}
\alias{BCI.env}
\docType{data}
\title{Barro Colorado Island Tree Counts}
\description{
  Tree counts in 1-hectare plots in the Barro Colorado Island and
  associated site information.
}
\usage{
data(BCI)
data(BCI.env)
}
\format{
  
  A data frame with 50 plots (rows) of 1 hectare with counts of trees on
  each plot with total of 225 species (columns). Full Latin names are
  used for tree species. The names were updated against
  \url{http://www.theplantlist.org} and Kress et al. (2009) which allows
  matching 207 of species against \doi{10.5061/dryad.63q27} (Zanne et
  al., 2014). The original species names are available as attribute
  \code{original.names} of \code{BCI}. See Examples for changed names.

  For \code{BCI.env}, a data frame with 50 plots (rows) and nine site
  variables derived from Pyke et al. (2001) and Harms et al. (2001):

  \describe{
    \item{\code{UTM.EW}: }{UTM coordinates (zone 17N) East-West.}
    \item{\code{UTM.NS}: }{UTM coordinates (zone 17N) North-South.}
    \item{\code{Precipitation}: }{Precipitation in mm per year.}
    \item{\code{Elevation}: }{Elevation in m above sea level.}
    \item{\code{Age.cat}: }{Forest age category.}
    \item{\code{Geology}: }{The Underlying geological formation.}

    \item{\code{Habitat}: }{Dominant habitat type based on the map of
    habitat types in 25 grid cells in each plot (Harms et al. 2001,
    excluding streamside habitat). The habitat types are \code{Young}
    forests (\emph{ca.} 100 years), old forests on > 7 degree slopes
    (\code{OldSlope}), old forests under 152 m elevation
    (\code{OldLow}) and at higher elevation (\code{OldHigh}) and
    \code{Swamp} forests.}

    \item{\code{River}: }{\code{"Yes"} if there is streamside habitat
    in the plot.}
    \item{\code{EnvHet}: }{Environmental Heterogeneity assessed as the
    Simpson diversity of frequencies of \code{Habitat} types in 25
    grid cells in the plot.}

}

 }
\details{

  Data give the numbers of trees at least 10 cm in diameter at breast
  height (DBH) in each one hectare quadrat in the 1982 BCI
  plot. Within each plot, all individuals were tallied and are
  recorded in this table. The full survey included smaller trees with
  DBH 1 cm or larger, but the \code{BCI} dataset is a subset of larger
  trees as compiled by Condit et al. (2002). The full data with
  thinner trees has densities above 4000 stems per hectare, or about
  ten times more stems than these data. The dataset \code{BCI} was
  provided (in 2003) to illustrate analysis methods in
  \pkg{vegan}. For scientific research on ecological issues we
  strongly recommend to access complete and more modern data (Condit
  et al. 2019) with updated taxonomy (Condit et al. 2020).

  The data frame contains only the Barro Colorado Island subset of the
  full data table of Condit et al. (2002).

  The quadrats are located in a regular grid. See \code{BCI.env} for the
  coordinates.

  A full description of the site information in \code{BCI.env} is
  given in Pyke et al. (2001) and Harms et al. (2001). \emph{N.B.}
  Pyke et al. (2001) and Harms et al. (2001) give conflicting
  information about forest age categories and elevation.

}
\source{

  \url{https://www.science.org/doi/10.1126/science.1066854} for
  community data and References for environmental data. For updated
  complete data (incl. thinner trees down to 1 cm), see Condit et
  al. (2019).

}

\seealso{

  Extra-CRAN package \pkg{natto}
  (\url{https://github.com/jarioksa/natto}) has data set
  \code{BCI.env2} with original grid data of Harms et al. (2001)
  habitat classification, and data set \code{BCI.taxon} of APG III
  classification of tree species.

}

\references{
  
  Condit, R, Pitman, N, Leigh, E.G., Chave, J., Terborgh, J., Foster,
  R.B., \enc{Nuñez}{Nunez}, P., Aguilar, S., Valencia, R., Villa, G.,
  Muller-Landau, H.C., Losos, E. & Hubbell, S.P. (2002).
  Beta-diversity in tropical forest trees. \emph{Science} 295,
  666--669.

  Condit R., Pérez, R., Aguilar, S., Lao, S., Foster, R. & Hubbell,
  S. (2019). Complete data from the Barro Colorado 50-ha plot: 423617
  trees, 35 years [Dataset].  \emph{Dryad}.
  \doi{https://doi.org/10.15146/5xcp-0d46}

  Condit, R., Aguilar, S., Lao, S., Foster, R., Hubbell,
  S. (2020). BCI 50-ha Plot Taxonomy [Dataset].  \emph{Dryad}.
  \doi{https://doi.org/10.15146/R3FH61}

  Harms K.E., Condit R., Hubbell S.P. & Foster R.B. (2001) Habitat
  associations of trees and shrubs in a 50-ha neotropical forest
  plot. \emph{J. Ecol.} 89, 947--959.

  Kress W.J., Erickson D.L, Jones F.A., Swenson N.G, Perez R., Sanjur
  O. & Bermingham E. (2009) Plant DNA barcodes and a community
  phylogeny of a tropical forest dynamics plot in Panama. \emph{PNAS}
  106, 18621--18626.

  Pyke, C. R., Condit, R., Aguilar, S., & Lao, S. (2001). Floristic
  composition across a climatic gradient in a neotropical lowland
  forest. \emph{Journal of Vegetation Science} 12, 553--566.
  \doi{10.2307/3237007}

  Zanne A.E., Tank D.C., Cornwell, W.K., Eastman J.M., Smith, S.A.,
  FitzJohn, R.G., McGlinn, D.J., O’Meara, B.C., Moles, A.T., Reich,
  P.B., Royer, D.L., Soltis, D.E., Stevens, P.F., Westoby, M., Wright,
  I.J., Aarssen, L., Bertin, R.I., Calaminus, A., Govaerts, R.,
  Hemmings, F., Leishman, M.R., Oleksyn, J., Soltis, P.S., Swenson,
  N.G., Warman, L. & Beaulieu, J.M. (2014) Three keys to the radiation
  of angiosperms into freezing environments. \emph{Nature} 506,
  89--92.  \doi{10.1038/nature12872} (published online Dec 22, 2013).

}

\examples{
data(BCI, BCI.env)
head(BCI.env)
## see changed species names
oldnames <- attr(BCI, "original.names")
taxa <- cbind("Old Names" = oldnames, "Current Names" = names(BCI))
noquote(taxa[taxa[,1] != taxa[,2], ])
}
\keyword{datasets}

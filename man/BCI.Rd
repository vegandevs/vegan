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
  
  A data frame with 50 plots (rows) of 1 hectare with counts of trees
  on each plot with total of 225 species (columns). Full Latin names
  are used for tree species. The names were updated against
  \url{http://www.theplantlist.org} in Jan 2014 (see ChangeLog 2.1-41
  for details) which allows matching 206 of species against
  \url{http://datadryad.org/resource/doi:10.5061/dryad.63q27} (Zanne
  et al., 2014).

  For \code{BCI.env}, a data frame with 50 plots (rows) and six site
  variables:

  \describe{
    \item{\code{UTM.EW}: }{UTM coordinates (zone 17N) East-West}
    \item{\code{UTM.NS}: }{UTM coordinates (zone 17N) North-South}
    \item{\code{Precipitation}: }{Precipitation in mm per year}
    \item{\code{Elevation}: }{Elevation in m above sea level}
    \item{\code{Age.cat}: }{Forest age category}
    \item{\code{Geology}: }{The Underlying geological formation}
  }

 }
\details{
  Data give the numbers of trees at least 10 cm in
  diameter at breast height (1.3 m above the ground) in each one hectare
  square of forest. Within each one hectare square, all individuals of
  all species were tallied and are recorded in this table.

  The data frame contains only the Barro Colorado Island subset of the
  original data.

  The quadrats are located in a regular grid. See \code{BCI.env} for the
  coordinates.

  A full description of the site information in \code{BCI.env} is given
  in Pyke et al (2001)

}
\source{
  \url{http://www.sciencemag.org/cgi/content/full/295/5555/666/DC1}
}

\references{
  
  Condit, R, Pitman, N, Leigh, E.G., Chave, J., Terborgh, J., Foster,
  R.B., \enc{Nuñez}{Nunez}, P., Aguilar, S., Valencia, R., Villa, G.,
  Muller-Landau, H.C., Losos, E. & Hubbell, S.P. (2002).
  Beta-diversity in tropical forest trees. \emph{Science} 295,
  666--669.

  Zanne A.E., Tank D.C., Cornwell, W.K., Eastman J.M., Smith, S.A.,
  FitzJohn, R.G., McGlinn, D.J., O’Meara, B.C., Moles, A.T., Reich,
  P.B., Royer, D.L., Soltis, D.E., Stevens, P.F., Westoby, M., Wright,
  I.J., Aarssen, L., Bertin, R.I., Calaminus, A., Govaerts, R.,
  Hemmings, F., Leishman, M.R., Oleksyn, J., Soltis, P.S., Swenson,
  N.G., Warman, L. & Beaulieu, J.M. (2014) Three keys to the radiation
  of angiosperms into freezing environments. \emph{Nature}
  doi:10.1038/nature12872 (published online Dec 22, 2013).

  Pyke, C. R., Condit, R., Aguilar, S., & Lao, S. (2001). Floristic
  composition across a climatic gradient in a neotropical lowland
  forest. \emph{Journal of Vegetation Science} 12, 553--566.
  doi:10.2307/3237007 

}
\examples{
data(BCI, BCI.env)
head(BCI.env)
}
\keyword{datasets}

\name{beetle}
\alias{beetle}
\alias{beetle.expl}
\docType{data}
\title{
Carabid beetle data with explanatory variables
}
\description{
Carabid beetle data collected at 192 sites by Colin Bergeron in the 
summer of 2003 in Northwestern Alberta, Canada. See Bergeron et al. 
(2011, 2012) and Blanchet et al. (2013) for details.
}
\usage{
data(beetle)
data(beetle.expl)
}
\format{
  The \code{beetle} data set is a data frame that contains the data on
  the following 37 carabid species:
  \describe{
    \item{\code{Agongrat}}{\emph{Agonum gratiosum}}
    \item{\code{Agonplac}}{\emph{Agonum placidum}}
    \item{\code{Agonretr}}{\emph{Agonum retractum}}
    \item{\code{Agonsord}}{\emph{Agonum sordens}}
    \item{\code{Agonsupe}}{\emph{Agonum superioris}}
    \item{\code{Amarlitt}}{\emph{Amara littoralis}}
    \item{\code{Amarluni}}{\emph{Amara lunicollis}}
    \item{\code{Badiobtu}}{\emph{Badister obtusus}}
    \item{\code{Bembgrap}}{\emph{Bembidion grapii}}
    \item{\code{Bembrupi}}{\emph{Bembidion rupicola}}
    \item{\code{Calaadve}}{\emph{Calathus advena}}
    \item{\code{Calaingr}}{\emph{Calathus ingratus}}
    \item{\code{Calofrig}}{\emph{Calosoma frigidum}}
    \item{\code{Caracham}}{\emph{Carabus chamissonis}}
    \item{\code{Dichcogn}}{\emph{Dicheirotrichus cognatus}}
    \item{\code{Elapamer}}{\emph{Elaphrus americanus}}
    \item{\code{Elaplapp}}{\emph{Elaphrus lapponicus}}
    \item{\code{Harpfulv}}{\emph{Harpalus fulvilabris}}
    \item{\code{Loripili}}{\emph{Loricera pilicornis}}
    \item{\code{Miscarct}}{\emph{Miscodera arctica}}
    \item{\code{Nebrgyll}}{\emph{Nebria gyllenhali}}
    \item{\code{Notibore}}{\emph{Notiophilus borealis}}
    \item{\code{Notidire}}{\emph{Notiophilus directus}}
    \item{\code{Patrfove}}{\emph{Patrobus foveocollis}}
    \item{\code{Patrsept}}{\emph{Patrobus septentrionis}}
    \item{\code{Platdece}}{\emph{Platynus decentis}}
    \item{\code{Platmann}}{\emph{Platynus mannerheimii}}
    \item{\code{Pteradst}}{\emph{Pterostichus adstrictus}}
    \item{\code{Pterbrev}}{\emph{Pterostichus brevicornis}}
    \item{\code{Pterpens}}{\emph{Pterostichus pensylvanicus}}
    \item{\code{Pterpunc}}{\emph{Pterostichus punctatissimus}}
    \item{\code{Pterripa}}{\emph{Pterostichus riparius}}
    \item{\code{Seriquad}}{\emph{Sericoda quadripunctata}}
    \item{\code{Sterhaem}}{\emph{Stereocerus haematopus}}
    \item{\code{Synuimpu}}{\emph{Synuchus impunctatus}}
    \item{\code{Trecapic}}{\emph{Trechus apicalis}}
    \item{\code{Trecchal}}{\emph{Trechus chalybeus}}

The \code{beetle.expl} data set is a data frame that contains the 
relative basal area of the 25 trees closest to the center of each 
sampled site. The relative basal area is presented by tree species. 

    \item{\code{Pt}}{Aspen (\emph{Populus tremuloides})}
    \item{\code{Bp}}{White birch (\emph{Betula papyrifera})}
    \item{\code{Ab}}{Balsam fir (\emph{Abie balsamea})}
    \item{\code{Ll}}{Tamarack (\emph{Larix laricina})}
    \item{\code{Pb}}{Balsam poplar (\emph{Populus balsamifera})}
    \item{\code{Pc}}{Lodgepole pine (\emph{Pinus contorta})}
    \item{\code{Pm}}{Black spruce (\emph{Picea mariana})}
    \item{\code{Pg}}{White spruce (\emph{Picea glauca})}
  }
}
\details{
In the \code{beetle} data, the abundance of carabids was divided by the 
number of days each trap was active in the field per sites. There were 
three traps per sites. This procedures was used to corrects for 
disturbances (such as destructions of traps), that may have occurred 
between visits.
}
\source{
All questions about these data should be adressed to Colin Bergeron 
(cb1[at]ualberta.ca).
}
\references{
Bergeron J.A.C., J.R. Spence, and W.J.A. Volney. 2011. Landscape 
patterns of species-level associations between ground-beetles 
(Coleoptera: Carabidae) and overstory trees in boreal forests of western
Canada (Coleoptera: Carabidae). In Erwin, TL (Ed), \emph{Proceedings of
a Symposium honoring the careers of Ross and Joyce Bell and their 
contributions to scientific work}, Burlington, VT, 12-15 June 2010. 
\emph{ZooKeys} \strong{147}: 577--600.

Bergeron J.A.C., F.G. Blanchet, J.R. Spence, and W.J.A. Volney. 2012. 
Ecosystem classification and inventory maps as surrogates for ground 
beetle assemblages in boreal forest. \emph{Journal of Plant Ecology} 
\strong{5}:97--108.

Blanchet, F.G., J.A.C. Bergeron, J.R. Spence, and F. He. 2013. Landscape
 effects of disturbance, habitat heterogeneity and spatial 
 autocorrelation for a ground beetle (Carabidae) assemblage in mature 
 boreal forest. \emph{Ecography} \strong{36}:636--647.

Blanchet, F.G., P. Legendre, J.A.C. Bergeron, and F. He. 2014. Consensus
RDA across dissimilarity coefficients for canonical ordination of 
community composition data. \emph{Ecological Monographs} 
\strong{84}: 491--511.
}
\examples{
data(beetle)
data(beetle.expl)
}
\keyword{datasets}

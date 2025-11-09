# Barro Colorado Island Tree Counts

Tree counts in 1-hectare plots in the Barro Colorado Island and
associated site information.

## Usage

``` r
data(BCI)
data(BCI.env)
```

## Format

A data frame with 50 plots (rows) of 1 hectare with counts of trees on
each plot with total of 225 species (columns). Full Latin names are used
for tree species. The names were updated with The Plant List web service
(now phased out) and Kress et al. (2009) which allows matching 207 of
species against
[doi:10.5061/dryad.63q27](https://doi.org/10.5061/dryad.63q27) (Zanne et
al., 2014). The original species names are available as attribute
`original.names` of `BCI`. See Examples for changed names.

For `BCI.env`, a data frame with 50 plots (rows) and nine site variables
derived from Pyke et al. (2001) and Harms et al. (2001):

- `UTM.EW`: :

  UTM coordinates (zone 17N) East-West.

- `UTM.NS`: :

  UTM coordinates (zone 17N) North-South.

- `Precipitation`: :

  Precipitation in mm per year.

- `Elevation`: :

  Elevation in m above sea level.

- `Age.cat`: :

  Forest age category.

- `Geology`: :

  The Underlying geological formation.

- `Habitat`: :

  Dominant habitat type based on the map of habitat types in 25 grid
  cells in each plot (Harms et al. 2001, excluding streamside habitat).
  The habitat types are `Young` forests (*ca.* 100 years), old forests
  on \> 7 degree slopes (`OldSlope`), old forests under 152 m elevation
  (`OldLow`) and at higher elevation (`OldHigh`) and `Swamp` forests.

- `River`: :

  `"Yes"` if there is streamside habitat in the plot.

- `EnvHet`: :

  Environmental Heterogeneity assessed as the Simpson diversity of
  frequencies of `Habitat` types in 25 grid cells in the plot.

## Details

Data give the numbers of trees at least 10 cm in diameter at breast
height (DBH) in each one hectare quadrat in the 1982 BCI plot. Within
each plot, all individuals were tallied and are recorded in this table.
The full survey included smaller trees with DBH 1 cm or larger, but the
`BCI` dataset is a subset of larger trees as compiled by Condit et al.
(2002). The full data with thinner trees has densities above 4000 stems
per hectare, or about ten times more stems than these data. The dataset
`BCI` was provided (in 2003) to illustrate analysis methods in vegan.
For scientific research on ecological issues we strongly recommend to
access complete and more modern data (Condit et al. 2019) with updated
taxonomy (Condit et al. 2020).

The data frame contains only the Barro Colorado Island subset of the
full data table of Condit et al. (2002).

The quadrats are located in a regular grid. See `BCI.env` for the
coordinates.

A full description of the site information in `BCI.env` is given in Pyke
et al. (2001) and Harms et al. (2001). *N.B.* Pyke et al. (2001) and
Harms et al. (2001) give conflicting information about forest age
categories and elevation.

## Source

<https://www.science.org/doi/10.1126/science.1066854> for community data
and References for environmental data. For updated complete data (incl.
thinner trees down to 1 cm), see Condit et al. (2019).

## See also

Extra-CRAN package natto (<https://github.com/jarioksa/natto>) has data
set `BCI.env2` with original grid data of Harms et al. (2001) habitat
classification, and data set `BCI.taxon` of APG III classification of
tree species.

## References

Condit, R, Pitman, N, Leigh, E.G., Chave, J., Terborgh, J., Foster,
R.B., Nuñez, P., Aguilar, S., Valencia, R., Villa, G., Muller-Landau,
H.C., Losos, E. & Hubbell, S.P. (2002). Beta-diversity in tropical
forest trees. *Science* 295, 666–669.

Condit R., Pérez, R., Aguilar, S., Lao, S., Foster, R. & Hubbell, S.
(2019). Complete data from the Barro Colorado 50-ha plot: 423617 trees,
35 years \[Dataset\]. *Dryad*.
[doi:10.15146/5xcp-0d46](https://doi.org/10.15146/5xcp-0d46)

Condit, R., Aguilar, S., Lao, S., Foster, R., Hubbell, S. (2020). BCI
50-ha Plot Taxonomy \[Dataset\]. *Dryad*.
[doi:10.15146/R3FH61](https://doi.org/10.15146/R3FH61)

Harms K.E., Condit R., Hubbell S.P. & Foster R.B. (2001) Habitat
associations of trees and shrubs in a 50-ha neotropical forest plot. *J.
Ecol.* 89, 947–959.

Kress W.J., Erickson D.L, Jones F.A., Swenson N.G, Perez R., Sanjur O. &
Bermingham E. (2009) Plant DNA barcodes and a community phylogeny of a
tropical forest dynamics plot in Panama. *PNAS* 106, 18621–18626.

Pyke, C. R., Condit, R., Aguilar, S., & Lao, S. (2001). Floristic
composition across a climatic gradient in a neotropical lowland forest.
*Journal of Vegetation Science* 12, 553–566.
[doi:10.2307/3237007](https://doi.org/10.2307/3237007)

Zanne A.E., Tank D.C., Cornwell, W.K., Eastman J.M., Smith, S.A.,
FitzJohn, R.G., McGlinn, D.J., O’Meara, B.C., Moles, A.T., Reich, P.B.,
Royer, D.L., Soltis, D.E., Stevens, P.F., Westoby, M., Wright, I.J.,
Aarssen, L., Bertin, R.I., Calaminus, A., Govaerts, R., Hemmings, F.,
Leishman, M.R., Oleksyn, J., Soltis, P.S., Swenson, N.G., Warman, L. &
Beaulieu, J.M. (2014) Three keys to the radiation of angiosperms into
freezing environments. *Nature* 506, 89–92.
[doi:10.1038/nature12872](https://doi.org/10.1038/nature12872)
(published online Dec 22, 2013).

## Examples

``` r
data(BCI, BCI.env)
head(BCI.env)
#>   UTM.EW  UTM.NS Precipitation Elevation Age.cat Geology  Habitat Stream EnvHet
#> 1 625754 1011569          2530       120      c3      Tb OldSlope    Yes 0.6272
#> 2 625754 1011669          2530       120      c3      Tb   OldLow    Yes 0.3936
#> 3 625754 1011769          2530       120      c3      Tb   OldLow     No 0.0000
#> 4 625754 1011869          2530       120      c3      Tb   OldLow     No 0.0000
#> 5 625754 1011969          2530       120      c3      Tb OldSlope     No 0.4608
#> 6 625854 1011569          2530       120      c3      Tb   OldLow     No 0.0768
## see changed species names
oldnames <- attr(BCI, "original.names")
taxa <- cbind("Old Names" = oldnames, "Current Names" = names(BCI))
noquote(taxa[taxa[,1] != taxa[,2], ])
#>       Old Names                     Current Names                 
#>  [1,] Abarema.macradenium           Abarema.macradenia            
#>  [2,] Acacia.melanoceras            Vachellia.melanoceras         
#>  [3,] Apeiba.aspera                 Apeiba.glabra                 
#>  [4,] Aspidosperma.cruenta          Aspidosperma.desmanthum       
#>  [5,] Cassipourea.elliptica         Cassipourea.guianensis        
#>  [6,] Cespedezia.macrophylla        Cespedesia.spathulata         
#>  [7,] Chlorophora.tinctoria         Maclura.tinctoria             
#>  [8,] Coccoloba.manzanillensis      Coccoloba.manzinellensis      
#>  [9,] Coussarea.curvigemmia         Coussarea.curvigemma          
#> [10,] Cupania.sylvatica             Cupania.seemannii             
#> [11,] Dipteryx.panamensis           Dipteryx.oleifera             
#> [12,] Eugenia.coloradensis          Eugenia.florida               
#> [13,] Eugenia.oerstedeana           Eugenia.oerstediana           
#> [14,] Guapira.standleyana           Guapira.myrtiflora            
#> [15,] Hyeronima.alchorneoides       Hieronyma.alchorneoides       
#> [16,] Inga.marginata                Inga.semialata                
#> [17,] Lonchocarpus.latifolius       Lonchocarpus.heptaphyllus     
#> [18,] Maquira.costaricana           Maquira.guianensis.costaricana
#> [19,] Phoebe.cinnamomifolia         Cinnamomum.triplinerve        
#> [20,] Swartzia.simplex.var.ochnacea Swartzia.simplex.continentalis
#> [21,] Tabebuia.guayacan             Handroanthus.guayacan         
```

# Oribatid Mite Data with Explanatory Variables

Oribatid mite data. 70 soil cores collected by Daniel Borcard in 1989.
See Borcard et al. (1992, 1994) for details.

## Usage

``` r
data(mite)
data(mite.env)
data(mite.pcnm)
data(mite.xy)
```

## Format

There are three linked data sets: `mite` that contains the data on 35
species of Oribatid mites, `mite.env` that contains environmental data
in the same sampling sites, `mite.xy` that contains geographic
coordinates, and `mite.pcnm` that contains 22 PCNM base functions
(columns) computed from the geographic coordinates of the 70 sampling
sites (Borcard & Legendre 2002). The whole sampling area was 2.5 m x 10
m in size.

The fields in the environmental data are:

- SubsDens:

  Substrate density (g/L)

- WatrCont:

  Water content of the substrate (g/L)

- Substrate:

  Substrate type, factor with levels
  `Sphagn1, Sphagn2 Sphagn3 Sphagn Litter Barepeat Interface`

- Shrub:

  Shrub density, an ordered factor with levels `1` \< `2` \< `3`

- Topo:

  Microtopography, a factor with levels `Blanket` and `Hummock`

## Source

Pierre Legendre

## References

Borcard, D., P. Legendre and P. Drapeau. 1992. Partialling out the
spatial component of ecological variation. Ecology 73: 1045-1055.

Borcard, D. and P. Legendre. 1994. Environmental control and spatial
structure in ecological communities: an example using Oribatid mites
(Acari, Oribatei). Environmental and Ecological Statistics 1: 37-61.

Borcard, D. and P. Legendre. 2002. All-scale spatial analysis of
ecological data by means of principal coordinates of neighbour matrices.
Ecological Modelling 153: 51-68.

## Examples

``` r
data(mite)
```

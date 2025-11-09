# Vegetation and Environment in Dutch Dune Meadows.

The dune meadow vegetation data, `dune`, has cover class values of 30
species on 20 sites. The corresponding environmental data frame
`dune.env` has following entries:

## Usage

``` r
data(dune)
  data(dune.env)
```

## Format

`dune` is a data frame of observations of 30 species at 20 sites. The
species names are abbreviated to 4+4 letters (see
[`make.cepnames`](https://vegandevs.github.io/vegan/reference/make.cepnames.md)).
The following names are changed from the original source (Jongman et al.
1987): *Leontodon autumnalis* to *Scorzoneroides*, and *Potentilla
palustris* to *Comarum*.

`dune.env` is a data frame of 20 observations on the following 5
variables:

- A1::

  a numeric vector of thickness of soil A1 horizon.

- Moisture::

  an ordered factor with levels: `1` \< `2` \< `4` \< `5`.

- Management::

  a factor with levels: `BF` (Biological farming), `HF` (Hobby farming),
  `NM` (Nature Conservation Management), and `SF` (Standard Farming).

- Use::

  an ordered factor of land-use with levels: `Hayfield` \< `Haypastu` \<
  `Pasture`.

- Manure::

  an ordered factor with levels: `0` \< `1` \< `2` \< `3` \< `4`.

## Source

Jongman, R.H.G, ter Braak, C.J.F & van Tongeren, O.F.R. (1987). *Data
Analysis in Community and Landscape Ecology*. Pudoc, Wageningen.

## Examples

``` r
data(dune)
data(dune.env)
```

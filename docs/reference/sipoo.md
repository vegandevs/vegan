# Birds in the Archipelago of Sipoo (Sibbo and Borgå)

Land birds on islands covered by coniferous forest in the Sipoo
Archipelago, southern Finland.

## Usage

``` r
data(sipoo)
  data(sipoo.map)
```

## Format

The `sipoo` data frame contains data of occurrences of 50 land bird
species on 18 islands in the Sipoo Archipelago (Simberloff & Martin,
1991, Appendix 3). The species are referred by 4+4 letter abbreviation
of their Latin names (but using five letters in two species names to
make these unique).

The `sipoo.map` data contains the geographic coordinates of the islands
in the ETRS89-TM35FIN coordinate system (EPSG:3067) and the areas of
islands in hectares.

## Source

Simberloff, D. & Martin, J.-L. (1991). Nestedness of insular avifaunas:
simple summary statistics masking complex species patterns. *Ornis
Fennica* 68:178–192.

## Examples

``` r
data(sipoo)
data(sipoo.map)
plot(N ~ E, data=sipoo.map, asp = 1)
```

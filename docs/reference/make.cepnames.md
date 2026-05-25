# Abbreviates a Two-Part Botanical or Zoological Latin Name into Character String

Function is based on
[`abbreviate`](https://rdrr.io/r/base/abbreviate.html), and will take
given number of characters from the first (genus) and last (epithet)
component of botanical or zoological Latin name and combine these into
one shorter character string. The names will be unique and more
characters will be used if needed. The default usage makes names with
4+4 characters popularized in Cornell Ecology Programs (CEP) and often
known as CEP names. Abbreviated names are useful in ordination plots and
other graphics to reduce clutter.

## Usage

``` r
make.cepnames(names, minlengths = c(4,4), seconditem = FALSE,
   uniqgenera = FALSE, named = FALSE, method)
```

## Arguments

- names:

  The names to be abbreviated into a vector abbreviated names.

- minlengths:

  The minimum lengths of first and second part of the abbreviation. If
  abbreviations are not unique, the parts can be longer.

- seconditem:

  Take always the second part of the original name to the abbreviated
  name instead of the last part.

- uniqgenera:

  Should the first part of the abbreviation (genus) also be unique.
  Unique genus can take space from the second part (epithet).

- method:

  The [`abbreviate`](https://rdrr.io/r/base/abbreviate.html) argument in
  last attempt to abbreviate the abbreviation. The default `method`
  tries to drop character from the end, but `"both.sides"` can remove
  characters from any position, including the genus part, and same genus
  can be abbreviated differently.

- named:

  Should the result vector be named by original `names`.

## Details

Cornell Ecology Programs (CEP) used eight-letter abbreviations for
species and site names. In species, the names were formed by taking four
first letters of the generic name and four first letters of the specific
or subspecific epithet. The current function produces CEP names as
default, but it can also use other lengths. The function is based on
[`abbreviate`](https://rdrr.io/r/base/abbreviate.html) and can produce
longer names if basic names are not unique. If generic name is shorter
than specified minimun length, more characters can be used by the
epithet. If `uniqgenera = TRUE` genus can use more characters, and these
reduce the number of characters available for the epithet. The function
drops characters from the end, but with `method = "both.sides"` the
function tries to drop characters from other positions, starting with
lower-case wovels, in the final attempt to abbreviate abbreviations.

## Value

Function returns a vector of abbreviated names.

## Author

Jari Oksanen

## See also

[`abbreviate`](https://rdrr.io/r/base/abbreviate.html).

## Note

The function does not handle Author names except strictly two-part names
with `seconditem = TRUE`. It is often useful to edit abbreviations
manually.

## Examples

``` r
names <- c("Aa maderoi", "Capsella bursa-pastoris", "Taraxacum",
  "Cladina rangiferina", "Cladonia rangiformis", "Cladonia cornuta",
  "Cladonia cornuta var. groenlandica", "Rumex acetosa",
  "Rumex acetosella")
make.cepnames(names)
#> [1] "Aamadero"    "Capsburs"    "Taraxacu"    "Cladrangife" "Cladrangifo"
#> [6] "Cladcorn"    "Cladgroe"    "Rumeacetosa" "Rumeacetose"
make.cepnames(names, uniqgenera = TRUE)
#> [1] "Aamadero"    "Capsburs"    "Taraxacu"    "Cladiran"    "Cladoran"   
#> [6] "Cladocor"    "Cladogro"    "Rumeacetosa" "Rumeacetose"
make.cepnames(names, method = "both.sides")
#> [1] "Aamadero" "Capsburs" "Taraxacu" "Cladrngf" "Cldrngfo" "Cladcorn" "Cladgroe"
#> [8] "Rumeacts" "Rmcetose"
```

# Reads a CEP (Canoco) data file

`read.cep` reads a file formatted with relaxed strict CEP format used in
Canoco software, among others.

## Usage

``` r
read.cep(file, positive=TRUE)
```

## Arguments

- file:

  File name (character variable).

- positive:

  Only positive entries, like in community data.

## Details

Cornell Ecology Programs (CEP) introduced several data formats designed
for punched cards. One of these was the ‘condensed strict’ format which
was adopted by popular software DECORANA and TWINSPAN. A relaxed variant
of this format was later adopted in Canoco software (ter Braak 1984).
Function `read.cep` reads legacy files written in this format.

The condensed CEP and CANOCO formats have:

- Two or three title cards, most importantly specifying the format and
  the number of items per record.

- Data in condensed format: First number on the line is the site
  identifier (an integer), and it is followed by pairs (‘couplets’) of
  numbers identifying the species and its abundance (an integer and a
  floating point number).

- Species and site names, given in Fortran format `(10A8)`: Ten names
  per line, eight columns for each.

With option `positive = TRUE` the function removes all rows and columns
with zero or negative marginal sums. In community data with only
positive entries, this removes empty sites and species. If data entries
can be negative, this ruins data, and such data sets should be read in
with option `positive = FALSE`.

## Value

Returns a data frame, where columns are species and rows are sites.
Column and row names are taken from the CEP file, and changed into
unique R names by [`make.names`](https://rdrr.io/r/base/make.names.html)
after stripping the blanks.

## References

ter Braak, C.J.F. (1984–): CANOCO – a FORTRAN program for *cano*nical
*c*ommunity *o*rdination by \[partial\] \[detrended\] \[canonical\]
correspondence analysis, principal components analysis and redundancy
analysis. *TNO Inst. of Applied Computer Sci., Stat. Dept. Wageningen,
The Netherlands*.

## Author

Jari Oksanen

## Note

Function `read.cep` used Fortran to read data in vegan 2.4-5 and
earlier, but Fortran I/O is no longer allowed in CRAN packages, and the
function was re-written in R. The original Fortran code was more robust,
and there are several legacy data sets that may fail with the current
version, but could be read with the previous Fortran version. CRAN
package cepreader makes available the original Fortran-based code run in
a separate subprocess. The cepreader package can also read ‘free’ and
‘open’ Canoco formats that are not handled in this function.

The function is based on
[`read.fortran`](https://rdrr.io/r/utils/read.fortran.html). If the
`REAL` format defines a decimal part for species abundances (such as
`F5.1`), [`read.fortran`](https://rdrr.io/r/utils/read.fortran.html)
divides the input with the corresponding power of 10 even when the input
data had explicit decimal separator. With `F5.1`, 100 would become 10,
and 0.1 become 0.01. Function `read.cep` tries to undo this division,
but you should check the scaling of results after reading the data, and
if necessary, multiply results to the original scale.

## Examples

``` r
## Provided that you have the file "dune.spe"
if (FALSE) { # \dontrun{
theclassic <- read.cep("dune.spe")} # }
```

# optspace: algorithm for matrix reconstruction from a partially revealed set

This function was adapted from the original source code in the Roptspace
R package (version 0.2.3; MIT License) by Raghunandan H. Keshavan,
Andrea Montanari, Sewoong Oh (2010). See `ROptSpace::OptSpace` for more
information. Let's assume an ideal matrix \\M\\ with \\(m\times n)\\
entries with rank \\r\\ and we are given a partially observed matrix
\\M\\E\\ which contains many missing entries. Matrix reconstruction - or
completion - is the task of filling in such entries. optspace is an
efficient algorithm that reconstructs \\M\\ from \\\|E\|=O(rn)\\
observed elements with relative root mean square error (RMSE) \$\$RMSE
\le C(\alpha)\sqrt{nr/\|E\|}\$\$.

## Usage

``` r
optspace(x, ropt = 3, niter = 5, tol = 1e-5, verbose = FALSE)
```

## Arguments

- x:

  An \\(n\times m)\\ matrix whose missing entries should be flagged as
  NA.

- ropt:

  `FALSE` to guess the rank, or a positive integer as a pre-defined rank
  (default: 3).

- niter:

  Maximum number of iterations allowed.

- tol:

  Stopping criterion for reconstruction in Frobenius norm.

- verbose:

  a logical value; `TRUE` to show progress, `FALSE` otherwise.

## Details

This implementation removes the trimming step of the original
`Roptspace::OptSpace` code in order to leave feature filtering to the
user. Some of the defaults have been adjusted to better reflect
ecological data. The implementation has been adjusted for ecological
applications as in Martino et al. (2019). The imputed matrix (\\M\\) in
the optspace output includes matrix reconstruction (\\XSY'\\).

## Value

Returns a named list containing:

- X:

  an \\(n \times r)\\ matrix as left singular vectors.

- S:

  an \\(r \times r)\\ matrix as singular values.

- Y:

  an \\(m \times r)\\ matrix as right singular vectors.

- dist:

  a vector containing reconstruction errors at each successive
  iteration.

- M:

  an \\(n \times m)\\ imputed matrix as calculated from `X`, `S` and
  `Y`.

## Author

Leo Lahti and Cameron Martino, with adaptations of the method
implemented in `Roptspace::OptSpace` by Keshavan et al. (2010).

## References

Keshavan, R. H., Montanari, A., Oh, S. (2010). Matrix Completion From a
Few Entries. *IEEE Transactions on Information Theory* **56**:2980–2998.
[doi:10.1109/TIT.2010.2046205](https://doi.org/10.1109/TIT.2010.2046205)
.

Martino, C., Morton, J.T., Marotz, C.A., Thompson, L.R., Tripathi, A.,
Knight, R. & Zengler, K. (2019) A novel sparse compositional technique
reveals microbial perturbations. *mSystems* **4**, 1.
[doi:10.1128/msystems.0016-19](https://doi.org/10.1128/msystems.0016-19)
.

## Examples

``` r
data(varespec)
# rclr transformation with no matrix completion for the 0/NA entries
x <- decostand(varespec, method = "rclr", impute = FALSE)
# Add matrix completion
xc <- optspace(x, ropt = 3, niter = 5, tol = 1e-5, verbose = FALSE)$M
```

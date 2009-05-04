\name{twostagechao}
\Rdversion{1.1}
\alias{twostagechao}
\title{
Multiple-Community Similarity Index
}
\description{
The function implements the two-stage probabilistic approach to 
multiple-community similarity indices proposed by Chao et al. (2008),
extending the traditional pairwise comparisons of sites.
}
\usage{
twostagechao(x, order = 2, N = nrow(x), m = 1, nboot = 200, subset)
}
\arguments{
  \item{x}{
A community matrix with rows as sites and columns as species.
}
  \item{order}{
The order of the similarity, see details.
}
  \item{N}{
The number of sites to use. By default it is the number of rows in \code{x}.
If \code{N < ncol(x)}, then \code{N} and \code{subset} must conform 
(\code{length(subset) == N}).
}
  \item{m}{
The exponent used for calculating the NESS index variant, see details.
}
  \item{nboot}{
Number of bootstrap samples to use for calculating standard error of the estimate.
If \code{nboot = 0}, standard error is not computed.
}
  \item{subset}{
Optional vector defining the subset of sites to use in the calculations.
If this argument is provided, its length should equal \code{N}.
}
}
\details{

PUT DETAILS HERE

}
\value{
The result is an object of class 'twostagechao' inheriting its structure from class 'htest'.
}
\references{
Chao, A., Jost, L., Chiang, S. C., Jiang, Y.-H., Chazdon, R. L. (2008)
A two-stage probabilistic approach to multiple-community similarity indices.
\emph{Biometrics} \bold{64}, 1178--1186.

Supplementary material at:
\url{http://www.biometrics.tibs.org/datasets/070539.pdf}
}
\author{
\enc{P\'eter S\'olymos}{Peter Solymos}, \email{solymos@ualberta.ca}
}
\seealso{
\code{\link{dist}}, \code{\link{vegdist}}
}
\examples{
data(BCI)
twostagechao(BCI)
res1 <- t(sapply(2:5, function(z) twostagechao(BCI, order=z)$statistic))
rownames(res1) <- paste("Order =", 2:5)
## similarity decreases as rare species get more weights
res1
## original example from Chao et al. 2008
## supplementary material, LSUR Primary Plot
lep1 <- structure(c(7L, 48L, 17L, 0L, 38L, 14L, 31L, 37L, 121L, 6L, 30L, 
16L, 0L, 27L, 10L, 0L, 24L, 5L, 2L, 23L, 6L, 3L, 21L, 6L, 1L, 
19L, 1L, 0L, 19L, 73L, 0L, 19L, 2L, 2L, 17L, 4L, 3L, 17L, 17L, 
2L, 16L, 5L, 0L, 16L, 11L, 0L, 16L, 4L, 2L, 15L, 1L, 0L, 15L, 
9L, 4L, 14L, 7L, 0L, 14L, 2L, 0L, 13L, 4L, 5L, 12L, 0L, 1L, 11L, 
7L, 2L, 11L, 20L, 1L, 9L, 0L, 2L, 9L, 5L, 0L, 9L, 0L, 1L, 8L, 
0L, 0L, 8L, 6L, 0L, 8L, 0L, 0L, 8L, 0L, 0L, 8L, 0L, 3L, 7L, 3L, 
0L, 7L, 17L, 0L, 7L, 0L, 1L, 0L, 0L, 1L, 0L, 0L, 1L, 0L, 0L, 
1L, 0L, 0L, 1L, 0L, 0L, 2L, 6L, 1L, 0L, 6L, 6L, 0L, 6L, 2L, 1L, 
5L, 0L, 2L, 5L, 1L, 1L, 5L, 3L, 0L, 5L, 2L, 0L, 5L, 1L, 0L, 5L, 
1L, 1L, 4L, 1L, 1L, 4L, 2L, 0L, 4L, 29L, 0L, 4L, 3L, 0L, 4L, 
1L, 0L, 4L, 0L, 0L, 3L, 2L, 0L, 3L, 1L, 0L, 3L, 0L, 3L, 2L, 0L, 
1L, 2L, 0L, 1L, 2L, 0L, 1L, 2L, 3L, 0L, 2L, 7L, 0L, 2L, 5L, 0L, 
2L, 1L, 0L, 2L, 1L, 0L, 2L, 1L, 0L, 2L, 0L, 0L, 2L, 0L, 0L, 2L, 
0L, 0L, 2L, 0L, 0L, 2L, 0L, 0L, 2L, 0L, 0L, 2L, 0L, 0L, 2L, 0L, 
3L, 0L, 1L, 1L, 0L, 1L, 1L, 0L, 3L, 0L, 0L, 54L, 0L, 0L, 3L, 
0L, 2L, 0L, 0L, 2L, 0L, 1L, 1L, 0L, 1L, 1L, 0L, 1L, 1L, 1L, 1L, 
1L, 2L, 0L, 1L, 5L, 0L, 1L, 3L, 0L, 1L, 3L, 0L, 1L, 1L, 0L, 1L, 
1L, 0L, 1L, 1L, 0L, 1L, 1L, 0L, 1L, 1L, 0L, 1L, 1L, 0L, 1L, 0L, 
0L, 1L, 0L, 0L, 1L, 0L, 0L, 1L, 0L, 0L, 1L, 0L, 0L, 1L, 0L, 0L, 
1L, 0L, 0L, 1L, 0L, 0L, 1L, 0L, 0L, 1L, 0L, 0L, 1L, 0L, 0L, 1L, 
0L, 0L, 1L, 0L, 0L, 1L, 0L, 0L, 1L, 0L, 0L, 1L, 0L, 0L, 1L, 0L, 
2L, 0L, 0L, 2L, 0L, 0L, 1L, 0L, 0L, 0L, 0L, 3L, 0L, 0L, 2L, 0L, 
0L, 1L, 0L, 0L, 1L, 0L, 0L, 1L), .Dim = c(3L, 120L), .Dimnames = list(
    c("Tree", "Sapling", "Seedling"), NULL))
## compare values in Table 3 of Chao et al. 2008
## (4th column, last 6 rows)
## results are different!
Order <- rep(2:3, 3)
m <- c(1, 1, 5, 5, 10, 10)
res2 <- t(sapply(1:6, function(z)
    twostagechao(lep1, order=Order[z], m=m[z])$statistic))
rownames(res2) <- c("Morisita C23", "Morisita C33", "NESS23(5)",
    "NESS33(5)", "NESS23(10)", "NESS33(10)")
res2
}
\keyword{ multivariate }

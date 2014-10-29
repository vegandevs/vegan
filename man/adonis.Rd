\encoding{UTF-8}
\name{adonis}
\alias{adonis}

\title{Permutational Multivariate Analysis of Variance Using Distance Matrices}

\description{Analysis of variance using distance matrices --- for
  partitioning distance matrices among sources of variation and fitting
  linear models (e.g., factors, polynomial regression) to distance 
  matrices; uses a permutation test with pseudo-\eqn{F} ratios.}

\usage{
adonis(formula, data, permutations = 999, method = "bray",
       strata = NULL, contr.unordered = "contr.sum",
       contr.ordered = "contr.poly", parallel = getOption("mc.cores"), ...)
}

\arguments{  
  \item{formula}{a typical model formula such as \code{Y ~ A + B*C}, but
    where \code{Y} is either a dissimilarity object (inheriting from
    class \code{"dist"}) or data frame or a matrix; \code{A}, \code{B}, and
    \code{C} may be factors or continuous variables. If a dissimilarity
  object is supplied, no species coefficients can be calculated (see
  Value below).} 
  \item{data}{ the data frame from which \code{A}, \code{B}, and
    \code{C} would be drawn.} 
  \item{permutations}{a list of control values for the permutations
    as returned by the function \code{\link[permute]{how}}, or the
    number of permutations required, or a permutation matrix where each
    row gives the permuted indices.}
  \item{method}{ the name of any method used in \code{\link{vegdist}} to
    calculate pairwise distances if the left hand side of the
    \code{formula} was a data frame or a matrix. } 
  \item{strata}{ groups (strata) within which to constrain permutations.  }
  \item{contr.unordered, contr.ordered}{contrasts used for the design
    matrix (default in R is dummy or treatment contrasts for unordered
    factors). }
  \item{parallel}{Number of parallel processes or a predefined socket
    cluster.  With \code{parallel = 1} uses ordinary, non-parallel
    processing. The parallel processing is done with \pkg{parallel}
    package.}
  \item{\dots}{Other arguments passed to \code{vegdist}.}
}

\details{\code{adonis} is a function for the analysis and partitioning
sums of squares using semimetric and metric distance matrices. Insofar
as it partitions sums of squares of a multivariate data set, it is
directly analogous to MANOVA (multivariate analysis of
variance). M.J. Anderson (McArdle and Anderson 2001, Anderson 2001) refers to the
method as \dQuote{permutational manova} (formerly \dQuote{nonparametric manova}). Further, as its inputs are
linear predictors, and a response matrix of an arbitrary number of
columns (2 to millions), it is a robust alternative to both parametric
MANOVA and to ordination methods for describing how variation is
attributed to different experimental treatments or uncontrolled
covariates. It is also analogous to redundancy analysis (Legendre and
Anderson 1999).

Typical uses of \code{adonis} include analysis of ecological community
data (samples X species matrices) or genetic data where we might have a
limited number of samples of individuals and thousands or millions of
columns of gene expression data (e.g. Zapala and Schork 2006).

\code{adonis} is an alternative to AMOVA (nested analysis of molecular
variance, Excoffier, Smouse, and Quattro, 1992;
\code{\link[ade4]{amova}} in the \pkg{ade4} package) for both crossed
and nested factors.

If the experimental design has nestedness, then use \code{strata} to
test hypotheses. For instance, imagine we are testing whether a
plant community is influenced by nitrate amendments, and we have two
replicate plots at each of two levels of nitrate (0, 10 ppm). We have
replicated the experiment in three fields with (perhaps) different
average productivity. In this design, we would need to specify
\code{strata = field} so that randomizations occur only \emph{within
each field} and not across all fields . See example below.

Like AMOVA (Excoffier et al. 1992), \code{adonis} relies on a
long-understood phenomenon that allows one to partition sums of squared
deviations from a centroid in two different ways (McArdle and Anderson
2001). The most widely recognized method, used, e.g., for ANOVA and
MANOVA, is to first identify the relevant centroids and then to
calculated the squared deviations from these points. For a centered
\eqn{n \times p}{n x p} response matrix \eqn{Y}, this method uses the
\eqn{p \times p}{p x p} inner product matrix \eqn{Y'Y}. The less
appreciated method is to use the \eqn{n \times n}{n x n} outer product
matrix \eqn{YY'}. Both AMOVA and \code{adonis} use this latter
method. This allows the use of any semimetric (e.g. Bray-Curtis, aka
Steinhaus, Czekanowski, and \enc{Sørensen}{Sorensen}) or metric
(e.g. Euclidean) distance matrix (McArdle and Anderson 2001). Using
Euclidean distances with the second method results in the same analysis
as the first method.

Significance tests are done using \eqn{F}-tests based on sequential sums
of squares from permutations of the raw data, and not permutations of
residuals. Permutations of the raw data may have better small sample
characteristics. Further, the precise meaning of hypothesis tests will
depend upon precisely what is permuted. The strata argument keeps groups
intact for a particular hypothesis test where one does not want to
permute the data among particular groups. For instance, \code{strata = B} 
causes permutations among levels of \code{A} but retains data within
levels of \code{B} (no permutation among levels of \code{B}). See
\code{\link{permutations}} for additional details on permutation tests
in Vegan.

The default \code{\link{contrasts}} are different than in \R in
general. Specifically, they use \dQuote{sum} contrasts, sometimes known
as \dQuote{ANOVA} contrasts. See a useful text (e.g. Crawley,
2002) for a transparent introduction to linear model contrasts. This
choice of contrasts is simply a personal
pedagogical preference. The particular contrasts can be set to any
\code{\link{contrasts}} specified in \R, including Helmert and treatment
contrasts.

Rules associated with formulae apply. See "An Introduction to R" for an
overview of rules.

\code{print.adonis} shows the \code{aov.tab} component of the output.
}
\value{
  This function returns typical, but limited, output for analysis of
  variance (general linear models). 

  \item{aov.tab}{Typical AOV table showing sources of variation,
    degrees of freedom, sequential sums of squares, mean squares,
    \eqn{F} statistics, partial \eqn{R^2}{R-squared} and \eqn{P}
    values, based on \eqn{N} permutations. }
  \item{coefficients}{ matrix of coefficients of the linear model, with
    rows representing sources of variation and columns representing
    species; each column represents a fit of a species abundance to the
    linear model. These are what you get when you fit one species to
    your predictors. These are NOT available if you supply the distance
    matrix in the formula, rather than the site x species matrix} 
  \item{coef.sites}{ matrix of coefficients of the linear model, with
    rows representing sources of variation and columns representing
    sites; each column represents a fit of a sites distances (from all
    other sites) to the  linear model. These are what you get when you
    fit distances of one site to
    your predictors. }   
  \item{f.perms}{ an \eqn{N} by \eqn{m} matrix of the null \eqn{F}
    statistics for each source of variation based on \eqn{N}
    permutations of the data. The permutations can be inspected with
    \code{\link{permustats}} and its support functions.}
  \item{model.matrix}{The \code{\link{model.matrix}} for the right hand
    side of the formula.}
  \item{terms}{The \code{\link{terms}} component of the model.}
}

\note{Anderson (2001, Fig. 4) warns that the method may confound
  location and dispersion effects: significant differences may be caused
  by different within-group variation (dispersion) instead of different
  mean values of the groups (see Warton et al. 2012 for a general
  analysis). However, it seems that \code{adonis} is less sensitive to
  dispersion effects than some of its alternatives (\code{\link{anosim}},
  \code{\link{mrpp}}). Function \code{\link{betadisper}} is a sister
  function to \code{adonis} to study the differences in dispersion
  within the same geometric framework.
}

\references{
Anderson, M.J. 2001. A new method for non-parametric multivariate
analysis of variance. \emph{Austral Ecology}, \strong{26}: 32--46.

Crawley, M.J. 2002. \emph{Statistical Computing: An Introduction to Data
  Analysis Using S-PLUS} 

Excoffier, L., P.E. Smouse, and J.M. Quattro. 1992. Analysis of
molecular variance inferred from metric distances among DNA haplotypes:
Application to human mitochondrial DNA restriction data. \emph{Genetics},
\strong{131}:479--491.

Legendre, P. and M.J. Anderson. 1999. Distance-based redundancy
analysis: Testing multispecies responses in multifactorial ecological
experiments. \emph{Ecological Monographs}, \strong{69}:1--24.

McArdle, B.H.  and M.J. Anderson. 2001. Fitting multivariate models to
community data: A comment on distance-based redundancy
analysis. \emph{Ecology}, \strong{82}: 290--297.

Warton, D.I., Wright, T.W., Wang, Y. 2012. Distance-based multivariate
analyses confound location and dispersion effects. \emph{Methods in
Ecology and Evolution}, 3, 89--101.

Zapala, M.A. and N.J. Schork. 2006. Multivariate regression analysis of
distance matrices for testing associations between gene expression
patterns and related variables. \emph{Proceedings of the National Academy of
Sciences, USA}, \strong{103}:19430--19435.
}
\author{Martin Henry H. Stevens
\email{HStevens@muohio.edu},
  adapted to \pkg{vegan} by Jari Oksanen. }

\seealso{ \code{\link{mrpp}}, \code{\link{anosim}},
  \code{\link{mantel}}, \code{\link{varpart}}. }
\examples{
data(dune)
data(dune.env)
adonis(dune ~ Management*A1, data=dune.env, permutations=99)


### Example of use with strata, for nested (e.g., block) designs.

dat <- expand.grid(rep=gl(2,1), NO3=factor(c(0,10)),field=gl(3,1) )
dat
Agropyron <- with(dat, as.numeric(field) + as.numeric(NO3)+2) +rnorm(12)/2
Schizachyrium <- with(dat, as.numeric(field) - as.numeric(NO3)+2) +rnorm(12)/2
total <- Agropyron + Schizachyrium
dotplot(total ~ NO3, dat, jitter.x=TRUE, groups=field,
        type=c('p','a'), xlab="NO3", auto.key=list(columns=3, lines=TRUE) )

Y <- data.frame(Agropyron, Schizachyrium)
mod <- metaMDS(Y)
plot(mod)
### Hulls show treatment
with(dat, ordihull(mod, group=NO3, show="0"))
with(dat, ordihull(mod, group=NO3, show="10", col=3))
### Spider shows fields
with(dat, ordispider(mod, group=field, lty=3, col="red"))

### Correct hypothesis test (with strata)
adonis(Y ~ NO3, data=dat, strata=dat$field, perm=999)

### Incorrect (no strata)
adonis(Y ~ NO3, data=dat, perm=999)
}

\keyword{multivariate }
\keyword{nonparametric }

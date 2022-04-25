\name{varpart}
\alias{varpart}
\alias{varpart2}
\alias{varpart3}
\alias{varpart4}
\alias{showvarparts}
\alias{summary.varpart}
\alias{plot.varpart}
\alias{plot.varpart234}
\alias{simpleRDA2}
\alias{simpleDBRDA}

\title{Partition the Variation of Community Matrix by 2, 3, or 4 Explanatory Matrices }

\description{ 

  The function partitions the variation in community data or community
  dissimilarities with respect to two, three, or four explanatory
  tables, using adjusted \eqn{R^2}{R-squared} in redundancy analysis
  ordination (RDA) or distance-based redundancy analysis. If response
  is a single vector, partitioning is by partial regression. Collinear
  variables in the explanatory tables do NOT have to be removed prior
  to partitioning.

}

\usage{
varpart(Y, X, ..., data, chisquare = FALSE, transfo, scale = FALSE,
    add = FALSE, sqrt.dist = FALSE, permutations)
\method{summary}{varpart}(object, ...)
showvarparts(parts, labels, bg = NULL, alpha = 63, Xnames,
    id.size = 1.2,  ...)
\method{plot}{varpart234}(x, cutoff = 0, digits = 1, ...)
}

\arguments{

\item{Y}{ Data frame or matrix containing the response data table or
  dissimilarity structure inheriting from \code{\link{dist}}. In
  community ecology, that table is often a site-by-species table or a
  dissimilarity object. }

\item{X}{Two to four explanatory models, variables or tables.  These can
  be defined in three alternative ways: (1) one-sided model formulae
  beginning with \code{~} and then defining the model, (2) name of a
  single numeric or factor variable, or (3) name of matrix with numeric
  or data frame with numeric and factor variables.  The model formulae
  can have factors, interaction terms and transformations of
  variables. The names of the variables in the model formula are found
  in data frame given in \code{data} argument, and if not found there,
  in the user environment.  Single variables, data frames or matrices
  are found in the user environment.  All entries till the next argument
  (\code{data} or \code{transfo}) are interpreted as explanatory models,
  and the names of these extra arguments cannot be abbreviated nor
  omitted.  }

\item{...}{Other parameters passed to functions. NB, arguments after
  dots cannot be abbreviated but they must be spelt out completely.}

\item{data}{The data frame with the variables used in the formulae in
  \code{X}.}

\item{chisquare}{Partition Chi-square or the inertia of Correspondence
  Analysis (\code{\link{cca}}).}

\item{transfo}{ Transformation for \code{Y} (community data) using
  \code{\link{decostand}}.  All alternatives in \code{decostand} can
  be used, and those preserving Euclidean metric include
  \code{"hellinger"}, \code{"chi.square"}, \code{"total"},
  \code{"norm"}. Ignored if \code{Y} are dissimilarities.}

\item{scale}{Should the columns of \code{Y} be standardized to unit
  variance. Ignored if \code{Y} are dissimilarities.}

\item{add}{Add a constant to the non-diagonal values to euclidify
  dissimilarities (see \code{\link{wcmdscale}} for details). Choice
  \code{"lingoes"} (or \code{TRUE}) use the recommended method of
  Legendre & Anderson (1999: \dQuote{method 1}) and \code{"cailliez"}
  uses their \dQuote{method 2}. The argument has an effect only when
  \code{Y} are dissimilarities.}

\item{sqrt.dist}{Take square root of dissimilarities. This often
  euclidifies dissimilarities. NB., the argument name cannot be
  abbreviated. The argument has an effect only when \code{Y} are
  dissimilarities.}

\item{permutations}{If \code{chisquare = TRUE}, the adjusted
    \eqn{R^2}{R-squared} is estimated by permutations, and this
    paramater can be a list of control values for the permutations as
    returned by the function \code{\link[permute]{how}}, or the number
    of permutations required, or a permutation matrix where each row
    gives the permuted indices.}

\item{parts}{Number of explanatory tables (circles) displayed.}
\item{labels}{Labels used for displayed fractions. Default is to use
  the same letters as in the printed output.}
\item{bg}{Fill colours of circles or ellipses.}
\item{alpha}{Transparency of the fill colour.  The argument takes
    precedence over possible transparency definitions of the
    colour. The value must be in range \eqn{0...255}, and low values
    are more transparent.  Transparency is not available in all
    graphics devices or file formats.}
  
\item{Xnames}{Names for sources of variation. Default names are \code{X1},
  \code{X2}, \code{X3} and \code{X4}. \code{Xnames=NA},
  \code{Xnames=NULL} and \code{Xnames=""} produce no names. The names
  can be changed to other names. It is often best to use short names. }

\item{id.size}{A numerical value giving the character expansion factor
  for the names of circles or ellipses. }

\item{x, object}{The \code{varpart} result.}
\item{cutoff}{The values below \code{cutoff} will not be displayed.}
\item{digits}{The number of significant digits; the number of decimal
  places is at least one higher.}
}

\details{

  The functions partition the variation in \code{Y} into components
  accounted for by two to four explanatory tables and their combined
  effects. If \code{Y} is a multicolumn data frame or matrix, the
  partitioning is based on redundancy analysis (RDA, see
  \code{\link{rda}}) or on constrained correspondence analysis if
  \code{chisquare = TRUE} (CCA, see \code{\link{cca}}).  If \code{Y}
  is a single variable, the partitioning is based on linear
  regression.  If \code{Y} are dissimilarities, the decomposition is
  based on distance-based redundancy analysis (db-RDA, see
  \code{\link{capscale}}) following McArdle & Anderson (2001). The
  input dissimilarities must be compatible to the results of
  \code{\link{dist}}. \pkg{Vegan} functions \code{\link{vegdist}},
  \code{\link{designdist}}, \code{\link{raupcrick}} and
  \code{\link{betadiver}} produce such objects, as do many other
  dissimilarity functions in \R packages. Partitioning will be made
  to squared dissimilarities analogously to using variance with
  rectangular data -- unless \code{sqrt.dist = TRUE} was specified.

  The function primarily uses adjusted \eqn{R^2}{R-squared} to assess
  the partitions explained by the explanatory tables and their
  combinations (see \code{\link{RsquareAdj}}), because this is the
  only unbiased method (Peres-Neto et al., 2006). The raw
  \eqn{R^2}{R-squared} for basic fractions are also displayed, but
  these are biased estimates of variation explained by the explanatory
  table. In correspondence analysis (\code{chisquare = TRUE}), the
  adjusted \eqn{R^2}{R-squared} are found by permutation and they vary
  in repeated analyses.

  The identifiable fractions are designated by lower case alphabets. The
  meaning of the symbols can be found in the separate document (use
  \code{browseVignettes("vegan")}), or can be displayed graphically
  using function \code{showvarparts}.

  A fraction is testable if it can be directly expressed as an RDA or
  db-RDA model.  In these cases the printed output also displays the
  corresponding RDA model using notation where explanatory tables
  after \code{|} are conditions (partialled out; see \code{\link{rda}}
  for details). Although single fractions can be testable, this does
  not mean that all fractions simultaneously can be tested, since the
  number of testable fractions is higher than the number of estimated
  models. The non-testable components are found as differences of
  testable components. The testable components have permutation
  variance in correspondence analysis (\code{chisquare = TRUE}), and
  the non-testable components have even higher variance.

  An abridged explanation of the alphabetic symbols for the individual
  fractions follows, but computational details should be checked in the
  vignette (readable with \code{browseVignettes("vegan")}) or in the
  source code.

  With two explanatory tables, the fractions explained 
  uniquely by each of the two tables are \code{[a]} and
  \code{[b]}, and their joint effect
  is  \code{[c]}.

  With three explanatory tables, the fractions explained uniquely
  by each of the three tables are  
  \code{[a]} to \code{[c]}, joint fractions between two tables are
  \code{[d]} to \code{[f]}, and the joint fraction between all three
  tables is \code{[g]}.

  With four explanatory tables, the fractions explained uniquely by each
  of the four tables are \code{[a]}
  to \code{[d]}, joint fractions between two tables are \code{[e]} to
  \code{[j]}, joint fractions between three variables are \code{[k]} to
  \code{[n]}, and the joint fraction between all four tables is
  \code{[o]}.

  \code{summary} will give an overview of unique and and overall
  contribution of each group of variables. The overall contribution
  (labelled as \dQuote{Contributed}) consists of the unique contribution
  of the variable and equal shares of each fraction where the variable
  contributes. The summary tabulates how each fraction is divided
  between the variables, and the contributed component is the sum of all
  these divided fractions. The summary is based on the idea of Lai et
  al. (2022), and is similar to the output of their \CRANpkg{rdacca.hp}
  package.

  There is a \code{plot} function that displays the Venn diagram and
  labels each intersection (individual fraction) with the adjusted R
  squared if this is higher than \code{cutoff}.  A helper function
  \code{showvarpart} displays the fraction labels. The circles and
  ellipses are labelled by short default names or by names defined by
  the user in argument \code{Xnames}. Longer explanatory file names can
  be written on the varpart output plot as follows: use option
  \code{Xnames=NA}, then add new names using the \code{text} function. A
  bit of fiddling with coordinates (see \code{\link{locator}}) and
  character size should allow users to place names of reasonably short
  lengths on the \code{varpart} plot.
  
}

\value{
  Function \code{varpart} returns an
  object of class \code{"varpart"} with items \code{scale} and
  \code{transfo} (can be missing) which hold information on
  standardizations, \code{tables} which contains names of explanatory
  tables, and \code{call} with the function \code{\link{call}}. The
  function \code{varpart} calls function \code{varpart2},
  \code{varpart3} or \code{varpart4} which return an object of class
  \code{"varpart234"} and saves its result in the item \code{part}.
  The items in this object are:
  \item{SS.Y}{Sum of squares of matrix \code{Y}.}
  \item{n}{Number of observations (rows).}
  \item{nsets}{Number of explanatory tables}
  \item{bigwarning}{Warnings on collinearity.}
  \item{fract}{Basic fractions from all estimated constrained models.}
  \item{indfract}{Individual fractions or all possible subsections in
    the Venn diagram (see \code{showvarparts}).}
  \item{contr1}{Fractions that can be found after conditioning on single
    explanatory table in models with three or four explanatory tables.}
  \item{contr2}{Fractions that can be found after conditioning on two
    explanatory tables in models with four explanatory tables.}
}
\section{Fraction Data Frames}{
  Items \code{fract},
  \code{indfract}, \code{contr1} and \code{contr2} are all data frames with
  items:
  \itemize{
  \item{\code{Df}: }{Degrees of freedom of numerator of the \eqn{F}-statistic
    for the fraction.}
  \item{\code{R.square}: }{Raw \eqn{R^2}{R-squared}. This is calculated only for
    \code{fract} and this is \code{NA} in other items.}
  \item{\code{Adj.R.square}: }{Adjusted \eqn{R^2}{R-squared}.}
  \item{\code{Testable}: }{If the fraction can be expressed as a (partial) RDA
    model, it is directly \code{Testable}, and this field is
    \code{TRUE}.  In that case the fraction label also gives the
    specification of the testable RDA model.}
  }
}

\references{ 

(a) References on variation partitioning

Borcard, D., P. Legendre & P. Drapeau. 1992. Partialling out the spatial
component of ecological variation. Ecology 73: 1045--1055.

Lai J., Y. Zou, J. Zhang & P. Peres-Neto. 2022. Generalizing
hierarchical and variation partitioning in multiple regression and
canonical analysis using the rdacca.hp R package. Methods in Ecology and
Evolution, 13: 782--788.

Legendre, P. & L. Legendre. 2012. Numerical ecology, 3rd English edition.
Elsevier Science BV, Amsterdam.

(b) Reference on transformations for species data

Legendre, P. and E. D. Gallagher. 2001. Ecologically meaningful
transformations for ordination of species data. Oecologia 129: 271--280.

(c) Reference on adjustment of the bimultivariate redundancy statistic

Peres-Neto, P., P. Legendre, S. Dray and D. Borcard. 2006. Variation partitioning
of species data matrices: estimation and comparison of fractions.
Ecology 87: 2614--2625.

(d) References on partitioning of dissimilarities

Legendre, P. & Anderson, M. J. (1999). Distance-based redundancy
analysis: testing multispecies responses in multifactorial ecological
experiments. \emph{Ecological Monographs} 69, 1--24.

McArdle, B.H. & Anderson, M.J. (2001). Fitting multivariate models
to community data: a comment on distance-based redundancy
analysis. Ecology 82, 290-297.

}

\author{ Pierre Legendre, Departement de Sciences Biologiques, Universite de
Montreal, Canada.  Further developed by Jari Oksanen. }

\note{

  You can use command \code{browseVignettes("vegan")} to display
  document which presents Venn diagrams showing the fraction names in
  partitioning the variation of Y with respect to 2, 3, and 4 tables of
  explanatory variables, as well as the equations used in variation
  partitioning.

  The functions frequently give negative estimates of variation.
  Adjusted \eqn{R^2}{R-squared} can be negative for any fraction;
  unadjusted \eqn{R^2}{R-squared} of testable fractions of variances
  will be non-negative.  Non-testable fractions cannot be found
  directly, but by subtracting different models, and these subtraction
  results can be negative.  The fractions are orthogonal, or linearly
  independent, but more complicated or nonlinear dependencies can
  cause negative non-testable fractions. Any fraction can be negative
  for non-Euclidean dissimilarities because the underlying db-RDA model
  can yield negative eigenvalues (see \code{\link{capscale}},
  \code{\link{dbrda}}). These negative eigenvalues in the underlying
  analysis can be avoided with arguments \code{sqrt.dist} and \code{add}
  which have a similar effect as in \code{\link{capscale}}: the square
  roots of several dissimilarities do not have negative eigenvalues, and
  no negative eigenvalues are produced after Lingoes or Cailliez
  adjustment, which in effect add random variation to the
  dissimilarities.

  A simplified, fast version of RDA, CCA adn dbRDA are used (functions
  \code{simpleRDA2}, \code{simpleCCA} and \code{simpleDBRDA}).  The
  actual calculations are done in functions \code{varpart2} to
  \code{varpart4}, but these are not intended to be called directly by
  the user.

}

\seealso{ 
  For analysing testable fractions, see \code{\link{rda}} and
  \code{\link{anova.cca}}. For data transformation, see
  \code{\link{decostand}}. Function \code{\link{inertcomp}} gives
  (unadjusted) components of variation for each species or site
  separately.  Function \code{\link{rda}} displays unadjusted
  components in its output, but \code{\link{RsquareAdj}} will give
  adjusted \eqn{R^2}{R-squared} that are similar to the current
  function also for partial models.  
}

\examples{
data(mite)
data(mite.env)
data(mite.pcnm)

# Two explanatory data frames -- Hellinger-transform Y
mod <- varpart(mite, mite.env, mite.pcnm, transfo="hel")
mod
summary(mod)

## Use fill colours
showvarparts(2, bg = c("hotpink","skyblue"))
plot(mod, bg = c("hotpink","skyblue"))
## Test fraction [a] using partial RDA, '~ .' in formula tells to use
## all variables of data mite.env.
aFrac <- rda(decostand(mite, "hel"), mite.env, mite.pcnm)
anova(aFrac)
## RsquareAdj gives the same result as component [a] of varpart
RsquareAdj(aFrac)

## Partition Bray-Curtis dissimilarities
varpart(vegdist(mite), mite.env, mite.pcnm)
## Three explanatory tables with formula interface
mod <- varpart(mite, ~ SubsDens + WatrCont, ~ Substrate + Shrub + Topo,
   mite.pcnm, data=mite.env, transfo="hel")
mod
summary(mod)
showvarparts(3, bg=2:4)
plot(mod, bg=2:4)

## Use RDA to test fraction [a]
## Matrix can be an argument in formula
rda.result <- rda(decostand(mite, "hell") ~ SubsDens + WatrCont +
   Condition(Substrate + Shrub + Topo) +
   Condition(as.matrix(mite.pcnm)), data = mite.env)
anova(rda.result)

## Four explanatory tables
mod <- varpart(mite, ~ SubsDens + WatrCont, ~Substrate + Shrub + Topo,
  mite.pcnm[,1:11], mite.pcnm[,12:22], data=mite.env, transfo="hel")
mod
summary(mod)
plot(mod, bg=2:5)
## Show values for all partitions by putting 'cutoff' low enough:
plot(mod, cutoff = -Inf, cex = 0.7, bg=2:5)
}

\keyword{ multivariate }

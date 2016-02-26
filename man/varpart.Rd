\name{varpart}
\alias{varpart}
\alias{varpart2}
\alias{varpart3}
\alias{varpart4}
\alias{showvarparts}
\alias{plot.varpart}
\alias{plot.varpart234}
\alias{simpleRDA2}

\title{Partition the Variation of Community Matrix by 2, 3, or 4 Explanatory Matrices }

\description{ 
  The function partitions the variation of response table Y with
  respect to two, three, or four explanatory tables, using adjusted
  \eqn{R^2}{R-squared} in redundancy analysis ordination (RDA). If Y
  contains a single vector, partitioning is by partial regression.
  Collinear variables in the explanatory tables do NOT have to be
  removed prior to partitioning.  
}

\usage{
varpart(Y, X, ..., data, transfo, scale = FALSE)
showvarparts(parts, labels, bg = NULL, alpha = 63, Xnames,
    id.size = 1.2,  ...)
\method{plot}{varpart234}(x, cutoff = 0, digits = 1, ...)
}

\arguments{
\item{Y}{ Data frame or matrix containing the response data
table. In community ecology, that table is often a site-by-species
table. }
\item{X}{Two to four explanatory models, variables or tables.  These can
  be defined in three alternative ways: (1) one-sided model formulae
  beginning with \code{~} and then defining the model, (2) name of a
  single numeric variable, or (3) name of data frame or matrix with
  numeric variables.  The model formulae can have factors,
  interaction terms and transformations of variables. The names of the
  variables in the model formula are found in data frame given in
  \code{data} argument, and if not found there, in the user
  environment.  Single numeric variables, data frames or matrices are
  found in the user environment.  All entries till the next argument
  (\code{data} or \code{transfo}) are interpreted as explanatory
  models, and the names of these arguments cannot be abbreviated nor
  omitted. 
  }
\item{data}{The data frame with the variables used in the formulae in
  \code{X}.} 
\item{transfo}{ Transformation for \code{Y} (community data) using
  \code{\link{decostand}}.  All alternatives in \code{decostand} can be
    used, and those preserving Euclidean metric include
    \code{"hellinger"}, \code{"chi.square"}, \code{"total"}, \code{"norm"}.}
\item{scale}{Should the columns of \code{Y} be standardized to unit variance}
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
  can be changed to other names. Use short names. }

\item{id.size}{A numerical value giving the character expansion factor
  for the names of circles or ellipses. }

\item{x}{The \code{varpart} result.}
\item{cutoff}{The values below \code{cutoff} will not be displayed.}
\item{digits}{The number of significant digits; the number of decimal
  places is at least one higher.}
\item{...}{Other parameters passed to functions.}
}

\details{
  The functions partition the variation in \code{Y} into components
  accounted for by two to four explanatory tables and their combined
  effects. If \code{Y} is a multicolumn data frame or
  matrix, the partitioning is based on redundancy analysis (RDA, see
  \code{\link{rda}}), and if \code{Y} is a single variable, the
  partitioning is based on linear regression.  

  The function primarily uses adjusted \eqn{R^2}{R-squared} to assess
  the partitions explained by the explanatory tables and their
  combinations, because this is the only unbiased method (Peres-Neto
  et al., 2006). The raw \eqn{R^2}{R-squared} for basic fractions are
  also displayed, but these are biased estimates of variation
  explained by the explanatory table. 

  The identifiable fractions are designated by lower case alphabets. The
  meaning of the symbols can be found in the separate document (use
  \code{browseVignettes("vegan")}), or can be displayed graphically
  using function \code{showvarparts}.

  A fraction is testable if it can be directly
  expressed as an RDA model.  In these cases the printed output also
  displays the corresponding RDA model using notation where explanatory
  tables after \code{|} are conditions (partialled out; see
  \code{\link{rda}} for details). Although single fractions can be
  testable, this does not mean that all fractions simultaneously can be
  tested, since there number of  testable fractions  is higher than
  the number of estimated models.

  An abridged explanation of the alphabetic symbols for the individual
  fractions follows, but computational details should be checked in the
  vignette (readable with \code{browseVignettes("vegan")}) or in the
  source code.

  With two explanatory tables, the fractions explained 
  uniquely by each of the two tables are \code{[a]} and
  \code{[c]}, and their joint effect
  is  \code{[b]} following Borcard et al. (1992).

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
  \item{Df}{Degrees of freedom of numerator of the \eqn{F}-statistic
    for the fraction.}
  \item{R.square}{Raw \eqn{R^2}{R-squared}. This is calculated only for
    \code{fract} and this is \code{NA} in other items.}
  \item{Adj.R.square}{Adjusted \eqn{R^2}{R-squared}.}
  \item{Testable}{If the fraction can be expressed as a (partial) RDA
    model, it is directly \code{Testable}, and this field is
    \code{TRUE}.  In that case the fraction label also gives the
    specification of the testable RDA model.}
  }
}

\references{ 

(a) References on variation partitioning

Borcard, D., P. Legendre & P. Drapeau. 1992. Partialling out the spatial
component of ecological variation. Ecology 73: 1045--1055.

Legendre, P. & L. Legendre. 2012. Numerical ecology, 3rd English edition.
Elsevier Science BV, Amsterdam.

(b) Reference on transformations for species data

Legendre, P. and E. D. Gallagher. 2001. Ecologically meaningful
transformations for ordination of species data. Oecologia 129: 271--280.

(c) Reference on adjustment of the bimultivariate redundancy statistic

Peres-Neto, P., P. Legendre, S. Dray and D. Borcard. 2006. Variation partitioning
of species data matrices: estimation and comparison of fractions.
Ecology 87: 2614--2625.
 }

\author{ Pierre Legendre, Departement de Sciences Biologiques, Universite de
Montreal, Canada.  Adapted to \pkg{vegan} by Jari Oksanen. }

\note{

  You can use command \code{browseVignettes("vegan")} to display
  document which presents Venn diagrams showing the fraction names in
  partitioning the variation of Y with respect to 2, 3, and 4 tables of
  explanatory variables, as well as the equations used in variation
  partitioning.

  The functions frequently give negative estimates of variation.
  Adjusted \eqn{R^2}{R-squared} can be negative for any fraction;
  unadjusted \eqn{R^2}{R-squared} of testable fractions always will be
  non-negative.  Non-testable fractions cannot be found directly, but
  by subtracting different models, and these subtraction results can
  be negative.  The fractions are orthogonal, or linearly independent,
  but more complicated or nonlinear dependencies can cause negative
  non-testable fractions.

  The current function will only use RDA in multivariate
  partitioning. It is much more complicated to estimate the adjusted
  R-squares for CCA, and unbiased analysis of CCA is not currently
  implemented.

  A simplified, fast version of RDA is used (function
  \code{simpleRDA2}).  The actual calculations are done in functions
  \code{varpart2} to \code{varpart4}, but these are not intended to be
  called directly by the user.

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

# Two explanatory matrices -- Hellinger-transform Y
# Formula shortcut "~ ." means: use all variables in 'data'.
mod <- varpart(mite, ~ ., mite.pcnm, data=mite.env, transfo="hel")
mod

## Use fill colours
showvarparts(2, bg = c("hotpink","skyblue"))
plot(mod, bg = c("hotpink","skyblue"))
# Alternative way of to conduct this partitioning
# Change the data frame with factors into numeric model matrix
mm <- model.matrix(~ SubsDens + WatrCont + Substrate + Shrub + Topo, mite.env)[,-1]
mod <- varpart(decostand(mite, "hel"), mm, mite.pcnm)
# Test fraction [a] using partial RDA:
aFrac <- rda(decostand(mite, "hel"), mm, mite.pcnm)
anova(aFrac, step=200, perm.max=200)
# RsquareAdj gives the same result as component [a] of varpart
RsquareAdj(aFrac)

# Three explanatory matrices 
mod <- varpart(mite, ~ SubsDens + WatrCont, ~ Substrate + Shrub + Topo,
   mite.pcnm, data=mite.env, transfo="hel")
mod
showvarparts(3, bg=2:4)
plot(mod, bg=2:4)
# An alternative formulation of the previous model using
# matrices mm1 amd mm2 and Hellinger transformed species data
mm1 <- model.matrix(~ SubsDens + WatrCont, mite.env)[,-1]
mm2 <- model.matrix(~ Substrate + Shrub + Topo, mite.env)[, -1]
mite.hel <- decostand(mite, "hel")
mod <- varpart(mite.hel, mm1, mm2, mite.pcnm)
# Use RDA to test fraction [a]
# Matrix can be an argument in formula
rda.result <- rda(mite.hel ~ mm1 + Condition(mm2) +
   Condition(as.matrix(mite.pcnm)))
anova(rda.result, step=200, perm.max=200)

# Four explanatory tables
mod <- varpart(mite, ~ SubsDens + WatrCont, ~Substrate + Shrub + Topo,
  mite.pcnm[,1:11], mite.pcnm[,12:22], data=mite.env, transfo="hel")
mod
plot(mod, bg=2:5)
# Show values for all partitions by putting 'cutoff' low enough:
plot(mod, cutoff = -Inf, cex = 0.7, bg=2:5)
}

\keyword{ multivariate }

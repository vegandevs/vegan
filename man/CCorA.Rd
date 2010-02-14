\name{CCorA}
\alias{CCorA}
\alias{print.CCorA}
\alias{biplot.CCorA}
\concept{ordination}

\title{Canonical Correlation Analysis}

\description{Canonical correlation analysis, following Brian McArdle's
unpublished graduate course notes, plus improvements to allow the
calculations in the case of very sparse and collinear matrices, and 
permutation test of Pillai's trace statistic. }

\usage{
CCorA(Y, X, stand.Y=FALSE, stand.X=FALSE, nperm = 0, ...)

\method{biplot}{CCorA}(x, plot.type="ov", xlabs, plot.axes = 1:2, int=0.5, 
   col.Y="red", col.X="blue", cex=c(0.7,0.9), ...)
}

\arguments{
  \item{Y}{ Left matrix (object class: \code{matrix} or \code{data.frame}). }
  \item{X}{ Right matrix (object class: \code{matrix} or \code{data.frame}). }
  \item{stand.Y}{ Logical; should \code{Y} be standardized? }
  \item{stand.X}{ Logical; should \code{X} be standardized? }
  \item{nperm}{ Numeric; number of permutations to evaluate the
    significance of Pillai's trace, e.g. \code{nperm=99} or \code{nperm=999}.}
  \item{x}{\code{CCoaR} result object.}
  \item{plot.type}{ A character string indicating which of the following 
    plots should be produced: \code{"objects"}, \code{"variables"}, \code{"ov"} 
    (separate graphs for objects and variables), or \code{"biplots"}. Any 
    unambiguous subset containing the first letters of these names can be used 
    instead of the full names. }
  \item{xlabs}{ Row labels. The default is to use row names, \code{NULL}
    uses row numbers instead, and \code{NA} suppresses plotting row names
    completely.}
  \item{plot.axes}{ A vector with 2 values containing the order numbers of 
    the canonical axes to be plotted. Default: first two axes. }
  \item{int}{ Radius of the inner circles plotted as visual references in 
    the plots of the variables. Default: \code{int=0.5}. With \code{int=0}, 
    no inner circle is plotted. }
  \item{col.Y}{ Color used for objects and variables in the first data 
    table (Y) plots. In biplots, the objects are in black. }
  \item{col.X}{ Color used for objects and variables in the second data 
    table (X) plots. }
  \item{cex}{ A vector with 2 values containing the size reduction factors 
    for the object and variable names, respectively, in the plots. 
    Default values: \code{cex=c(0.7,0.9)}. }
  \item{\dots}{ Other arguments passed to these functions. The function 
    \code{biplot.CCorA} passes graphical arguments to \code{\link{biplot}} 
    and \code{\link{biplot.default}}. \code{CCorA} currently ignores extra 
    arguments. } 
}

\details{
  Canonical correlation analysis (Hotelling 1936) seeks linear
  combinations of the variables of \code{Y} that are maximally
  correlated to linear combinations of the variables of \code{X}. The
  analysis estimates the relationships and displays them in graphs.
  Pillai's trace statistic is computed and tested parametrically (F-test);
  a permutation test is also available.

  Algorithmic note -- 
  The blunt approach would be to read the two matrices, compute the
  covariance matrices, then the matrix
  \code{S12 \%*\% inv(S22) \%*\% t(S12) \%*\% inv(S11)}.
  Its trace is Pillai's trace statistic. 
  This approach may fail, however, when there is heavy multicollinearity
  in very sparse data matrices. The safe approach is to replace all data
  matrices by their PCA object scores.

The function can produce different types of plots depending on the option 
chosen: 
\code{"objects"} produces two plots of the objects, one in the space 
of Y, the second in the space of X; 
\code{"variables"} produces two plots of the variables, one of the variables 
of Y in the space of Y, the second of the variables of X in the space of X; 
\code{"ov"} produces four plots, two of the objects and two of the variables; 
\code{"biplots"} produces two biplots, one for the first matrix (Y) and 
one for second matrix (X) solutions. For biplots, the function passes all arguments 
to \code{\link{biplot.default}}; consult its help page for configuring biplots.
}

\value{
  Function \code{CCorA} returns a list containing the following elements:
  \item{ Pillai }{ Pillai's trace statistic = sum of the canonical
    eigenvalues. } 
  \item{ Eigenvalues }{ Canonical eigenvalues. They are the squares of the
    canonical correlations. }
  \item{ CanCorr }{ Canonical correlations. }
  \item{ Mat.ranks }{ Ranks of matrices \code{Y} and \code{X}. }
  \item{ RDA.Rsquares }{ Bimultivariate redundancy coefficients
    (R-squares) of RDAs of Y|X and X|Y. } 
  \item{ RDA.adj.Rsq }{ \code{RDA.Rsquares} adjusted for \code{n} and the number 
    of explanatory variables. }
  \item{ nperm }{ Number of permutations. }
  \item{ p.Pillai }{ Parametric probability value associated with Pillai's trace. }
  \item{ p.perm }{ Permutational probability associated with Pillai's trace. }
  \item{ Cy }{ Object scores in Y biplot. }
  \item{ Cx }{ Object scores in X biplot. }
  \item{ corr.Y.Cy }{ Scores of Y variables in Y biplot, computed as cor(Y,Cy). }
  \item{ corr.X.Cx }{ Scores of X variables in X biplot, computed as cor(X,Cx). }
  \item{ corr.Y.Cx }{ cor(Y,Cy) available for plotting variables Y in space of X manually. }
  \item{ corr.X.Cy }{ cor(X,Cx) available for plotting variables X in space of Y manually. }
  \item{ call }{ Call to the CCorA function. }
}

\references{ 
  Hotelling, H. 1936. Relations between two sets of
  variates. \emph{Biometrika} \strong{28}: 321-377.
  
  Legendre, P. 2005. Species associations: the Kendall coefficient of 
  concordance revisited. \emph{Journal of Agricultural, Biological, and 
  Environmental Statistics} \strong{10}: 226-245.
}

\author{ Pierre Legendre, Departement de Sciences Biologiques,
  Universite de Montreal. Implemented in \pkg{vegan} with the help of
  Jari Oksanen. }

\examples{
# Example using two mite groups. The mite data are available in vegan
data(mite)
# Two mite species associations (Legendre 2005, Fig. 4)
group.1 <- c(1,2,4:8,10:15,17,19:22,24,26:30)
group.2 <- c(3,9,16,18,23,25,31:35)
# Separate Hellinger transformations of the two groups of species 
mite.hel.1 <- decostand(mite[,group.1], "hel")
mite.hel.2 <- decostand(mite[,group.2], "hel")
rownames(mite.hel.1) = paste("S",1:nrow(mite),sep="")
rownames(mite.hel.2) = paste("S",1:nrow(mite),sep="")
out <- CCorA(mite.hel.1, mite.hel.2)
out
summary(out)
biplot(out, "ob")                 # Two plots of objects
biplot(out, "v", cex=c(0.7,0.6))  # Two plots of variables
biplot(out, "ov", cex=c(0.7,0.6)) # Four plots (2 for objects, 2 for variables)
biplot(out, "b", cex=c(0.7,0.6))  # Two biplots
biplot(out, xlabs = NA, plot.axes = c(3,5))    # Plot axes 3, 5. No object names
biplot(out, plot.type="biplots", xlabs = NULL) # Replace object names by numbers

# Example using random numbers. No significant relationship is expected
mat1 <- matrix(rnorm(60),20,3)
mat2 <- matrix(rnorm(100),20,5)
out2 = CCorA(mat1, mat2, nperm=99)
out2
biplot(out2, "b")
}

\keyword{ multivariate }

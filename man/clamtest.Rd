\name{clamtest}
\alias{clamtest}
\alias{summary.clamtest}
\alias{print.summary.clamtest}
\alias{plot.clamtest}
\title{
Multinomial Species Classification Method (CLAM)
}
\description{
The CLAM statistical approach for classifying generalists and
specialists in two distinct habitats is described in Chazdon et al. (2011).
}
\usage{
clamtest(comm, groups, coverage.limit = 10, specialization = 0.667, 
npoints = 20, alpha = 0.05/npoints, ...)
\method{summary}{clamtest}(object, ...)
\method{print}{summary.clamtest}(x, 
digits = max(3, getOption("digits") - 3), ...)
\method{plot}{clamtest}(x, xlab, ylab, main, 
pch = 21:24, col.points = 1:4, col.lines = 2:4, 
lty = 1:3, position = "bottomright", ...)
}
\arguments{
  \item{comm}{
Community matrix, consisting of counts.
}
  \item{groups}{
A vector identifying the two habitats. Must have exactly
two unique values or levels.
}
  \item{coverage.limit}{
Integer, below this limit the sample coverage based correction
is applied to rare species. Sample coverage is claculated separately 
for the two habitats. Sample relative abundances are used for species 
with higher than or equal to 10 total counts per habitat.
}
  \item{specialization}{
Numeric, specialization threshold value between 0 and 1.
The value of 0.667 represents 'supermajority' rule,
while a value of 0.5 represents a 'simple majority' rule
to assign shared species as habitat specialists.
}
  \item{npoints}{
Integer, number of points used to determine the boundary lines
in the plots.
}
  \item{alpha}{
Numeric, significance level for individual tests.
The default value sets the overall (experiment-wise) 
significance level to 0.05 by controlling for
the number of points (\code{npoints}) used to construct
boundary lines.
}
  \item{x, object}{
Fitted model object of class \code{"clamtest"}.
}
  \item{xlab, ylab}{
Labels for the plot axes.
}
  \item{main}{
Main title of the plot.
}
  \item{pch, col.points}{
Symbols and color used in plotting species groups.
}
  \item{lty, col.lines}{
Line type and color for boundary lines in plot to separate species groups.
}
  \item{position}{
Position of figure legend, see \code{\link{legend}} for description.
Legend not shown if \code{position = NULL}.
}
  \item{\dots}{
Additional arguments passed to methods.
}
}
\details{
The method uses a multinomial model based on estimated species
relative abundance in two habitats (A, B), it minimizes bias due to differences in sampling
intensities between two habitat types as well as bias due to insufficient sampling within each
habitat. The method permits a robust statistical classification of habitat specialists and
generalists, without excluding rare species a priori (Chazdon et al. 2011).
Based on a user-defined \code{specialization}
threshold, the model classifies species into one of four groups: (1) generalist; (2) habitat A
specialist; (3) habitat B specialist; and (4) too rare to classify with confidence.
}
\value{
A data frame (with class attribute \code{"clamtest"}), with columns:
\itemize{
  \item{\code{Species}:}{ species names (column names from \code{comm}),}
  \item{\code{Total_*A*}:}{ total count in habitat A,}
  \item{\code{Total_*B*}:}{ total count in habitat B,}
  \item{\code{Classes}:}{ species classification, a factor with levels
    \code{Generalist}, \code{Specialist_*A*}, \code{Specialist_*B*}, and \code{Too_rare}.}
}
\code{*A*} and \code{*B*} are placeholders for habitat names/labels found in the data.

The \code{summary} method returns descriptive statistics of the results.
The \code{plot} method returns values invisibly and produces a bivariate
scatterplot of species total abundances in the two habitats. Symbols and
boundary lines are shown to separate species groups.
}
\references{
Chazdon, R. L., Chao, A., Colwell, R. K., Lin, S.-Y., Norden, N., 
Letcher, S. G., Clark, D. B., Finegan, B. and Arroyo J. P.(2011). 
A novel statistical method for classifying habitat
generalists and specialists. \emph{Ecology} \bold{92}, 1332--1343.
}
\author{
Peter Solymos \email{solymos@ualberta.ca}
}
\note{
The code was tested against standalone CLAM software provided
on the website of Anne Chao (\url{http://chao.stat.nthu.edu.tw/softwarece.html});
Minor inconsistencies were found, especially for finding the
threshold for 'too rare' species.
These inconsistencies are probably due to numerical differences between the
two implementation. The current R implementation uses 
root finding for iso-lines instead of itarative search.
}
\examples{
data(mite)
data(mite.env)
x <- clamtest(mite, mite.env$Shrub=="None", alpha=0.005, specialization = 0.667)
summary(x)
head(x)
plot(x)
}
\keyword{ htest }

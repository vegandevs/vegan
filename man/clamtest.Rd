\name{clamtest}
\alias{clamtest}
\alias{summary.clamtest}
\alias{plot.clamtest}

\title{
Multinomial Species Classification Method (CLAM)
}

\description{
The CLAM statistical approach for classifying generalists and
specialists in two distinct habitats is described in Chazdon et al. (2011).
}
\usage{
clamtest(comm, groups, coverage.limit = 10, specialization = 2/3, 
   npoints = 20, alpha = 0.05/20)
\method{summary}{clamtest}(object, ...)
\method{plot}{clamtest}(x, xlab, ylab, main,  pch = 21:24, col.points = 1:4, 
   col.lines = 2:4, lty = 1:3, position = "bottomright", ...)
}
\arguments{
  \item{comm}{
Community matrix, consisting of counts.
}
  \item{groups}{
A vector identifying the two habitats. Must have exactly
two unique values or levels. Habitat IDs in the grouping vector
must match corresponding rows in the community matrix \code{comm}.
}
  \item{coverage.limit}{
Integer, the sample coverage based correction 
is applied to rare species with counts below this limit. 
Sample coverage is calculated separately 
for the two habitats. Sample relative abundances are used for species 
with higher than or equal to \code{coverage.limit} total counts per habitat.
}
  \item{specialization}{
Numeric, specialization threshold value between 0 and 1.
The value of \eqn{2/3} represents \sQuote{supermajority} rule,
while a value of \eqn{1/2} represents a \sQuote{simple majority} rule
to assign shared species as habitat specialists.
}
  \item{npoints}{
Integer, number of points used to determine the boundary lines
in the plots.
}
\item{alpha}{ Numeric, nominal significance level for individual
  tests.  The default value reduces the conventional limit of
  \eqn{0.05} to account for overdispersion and multiple testing for
  several species simultaneously. However, the is no firm reason for
  exactly this limit.  }
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
Symbols and colors used in plotting species groups.
}
  \item{lty, col.lines}{
Line types and colors for boundary lines in plot to separate species groups.
}
  \item{position}{
Position of figure legend, see \code{\link{legend}} for specification details.
Legend not shown if \code{position = NULL}.
}
  \item{\dots}{
Additional arguments passed to methods.
}
}

\details{ The method uses a multinomial model based on estimated
  species relative abundance in two habitats (A, B). It minimizes bias
  due to differences in sampling intensities between two habitat types
  as well as bias due to insufficient sampling within each
  habitat. The method permits a robust statistical classification of
  habitat specialists and generalists, without excluding rare species
  \emph{a priori} (Chazdon et al. 2011).  Based on a user-defined
  \code{specialization} threshold, the model classifies species into
  one of four groups: (1) generalists; (2) habitat A specialists; (3)
  habitat B specialists; and (4) too rare to classify with confidence.
  } 

\value{ A data frame (with class attribute \code{"clamtest"}),
  with columns: 

    \item{\code{Species}:}{ species name (column names from \code{comm}),} 
    \item{\code{Total_*A*}:}{ total count in habitat A,} 
    \item{\code{Total_*B*}:}{ total count in habitat B,} 
    \item{\code{Classes}:}{ species classification, a factor with
       levels \code{Generalist}, \code{Specialist_*A*},
       \code{Specialist_*B*}, and \code{Too_rare}.}  

  \code{*A*} and \code{*B*} are placeholders for habitat names/labels found in the
  data.

The \code{summary} method returns descriptive statistics of the results.
The \code{plot} method returns values invisibly and produces a bivariate
scatterplot of species total abundances in the two habitats. Symbols and
boundary lines are shown for species groups.
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
on the website of Anne Chao (which were then at http://chao.stat.nthu.edu.tw/wordpress);
minor inconsistencies were found, especially for finding the
threshold for 'too rare' species.
These inconsistencies are probably due to numerical differences between the
two implementation. The current \R implementation uses 
root finding for iso-lines instead of iterative search.

The original method (Chazdon et al. 2011) has two major problems:
\enumerate{
  
  \item It assumes that the error distribution is multinomial. This is
    a justified choice if individuals are freely distributed, and
    there is no over-dispersion or clustering of individuals. In most
    ecological data, the variance is much higher than multinomial
    assumption, and therefore test statistic are too optimistic.

  \item The original authors suggest that multiple testing adjustment
    for multiple testing should be based on the number of points
    (\code{npoints}) used to draw the critical lines on the plot,
    whereas the adjustment should be based on the number of tests (i.e.,
    tested species). The function uses the same numerical values as
    the original paper, but there is no automatic connection between
    \code{npoints} and \code{alpha} arguments, but you must work out
    the adjustment yourself.
}
}
\examples{
data(mite)
data(mite.env)
sol <- with(mite.env, clamtest(mite, Shrub=="None", alpha=0.005))
summary(sol)
head(sol)
plot(sol)
}
\keyword{ htest }

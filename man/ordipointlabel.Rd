\name{ordipointlabel}
\alias{ordipointlabel}
\alias{plot.ordipointlabel}

\title{ Ordination Plots with Points and Optimized Locations for Text }
\description{
  The function \code{ordipointlabel} produces ordination plots with
  points and text label to the points. The points are in the exact
  location given by the ordination, but the function tries to optimize
  the location of the text labels to minimize overplotting text. The
  function may be useful with moderately crowded ordination plots.
}
\usage{
ordipointlabel(x, display = c("sites", "species"), choices = c(1, 2),
   col = c(1, 2),  pch = c("o", "+"), font = c(1, 1), 
   cex = c(0.8, 0.8), add = FALSE, select, ...)

\method{plot}{ordipointlabel}(x, ...)
}

\arguments{
  \item{x}{For \code{ordipointlabel()} a result object from an
    ordination function. For \code{plot.ordipointlabel} an object
    resulting from a call to \code{ordipointlabel()}.}
  \item{display}{Scores displayed in the plot. }
  \item{choices}{Axes shown. }
  \item{col, pch, font, cex}{Colours, point types, font style and
     character expansion for each kind of scores displayed in the
     plot. These should be vectors of the same length as the number of
     items in \code{display}.}
  \item{add}{ Add to an existing plot. }
  \item{select}{Items to be displayed.  This can either be a logical
    vector which is \code{TRUE} for displayed items or a vector of indices
    of displayed items. \code{select} is only used if a single set of
    scores is being plotted (i.e. \code{length(display) == 1}),
    otherwise it is ignored and a warning issued. If a logical vector is
    used, it must have the same length as the scores plotted.}
  \item{\dots}{Other arguments passed to \code{\link{points}} and
  \code{\link{text}}.}
}
\details{
  The function uses simulated annealing (\code{\link{optim}},
  \code{method = "SANN"}) to optimize the location of the text labels
  to the points. There are eight possible locations: up, down, sides
  and corners. There is a weak preference to text right above the
  point, and a weak avoidance of corner positions. The exact locations
  and the goodness of solution varies between runs, and there is no
  guarantee of finding the global optimum. The optimization can take a
  long time in difficult cases with a high number of potential
  overlaps. Several sets of scores can be displayed in one plot. 

  The function is modelled after \code{\link[maptools]{pointLabel}} in
  \pkg{maptools} package (which has chained dependencies of \code{S4}
  packages). 
}
\value{
  The function returns invisibly an object of class
  \code{ordipointlabel} with items \code{xy} for coordinates of
  points, \code{labels} for coordinates of labels, items \code{pch},
  \code{cex} and \code{font} for graphical parameters of each point or
  label. In addition, it returns the result of \code{\link{optim}} as
  an attribute \code{"optim"}. The unit of overlap is the area
  of character \code{"m"}, and with variable \code{cex} it is the
  smallest alternative.

  There is a \code{plot} method based on \code{orditkplot} but which
  does not alter nor reset the graphical parameters via \code{par}.
  
  The result object from \code{ordipointlabel} inherits from
  \code{\link{orditkplot}}, and can also be replotted with its
  \code{plot} method. It may be possible to further edit the result
  object with \code{\link{orditkplot}}, but for good results it is
  necessary that the points span the whole horizontal axis without empty
  margins. 
}
\references{ See \code{\link[maptools]{pointLabel}} for references. }
\author{ Jari Oksanen }
\note{ 
  The function is designed for ordination graphics, and the
  optimization works properly with plots of isometric aspect ratio.
}
\seealso{ \code{\link[maptools]{pointLabel}} for the model
  implementation, and \code{\link{optim}} for the optimization. }
\examples{
data(dune)
ord <- cca(dune)
plt <- ordipointlabel(ord)

## set scaling - should be no warnings!
ordipointlabel(ord, scaling = 1)

## plot then add
plot(ord, scaling = 3, type = "n")
ordipointlabel(ord, display = "species", scaling = 3, add = TRUE)
ordipointlabel(ord, display = "sites", scaling = 3, add = TRUE)

## redraw plot without rerunning SANN optimisation
plot(plt)
}
\keyword{ hplot }
\keyword{ aplot }

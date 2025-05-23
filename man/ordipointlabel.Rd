\name{ordipointlabel}
\alias{ordipointlabel}
\alias{plot.ordipointlabel}
\alias{labels.ordipointlabel}

\title{ Ordination Plots with Points and Optimized Locations for Text }

\description{ Function produces ordination plots with points and text
  labels to the points. The points are in the fixed locations given by
  the ordination, but the locations of the text labels are optimized to
  minimize overplotting. The function is useful with moderately crowded
  ordination plots.
}

\usage{
ordipointlabel(x, display = c("sites", "species"), choices = c(1, 2),
   col = c(1, 2),  pch = c("o", "+"), font = c(1, 1), 
   cex = c(0.7, 0.7), add = inherits(x, "ordiplot"), labels, bg, select, ...)

\method{plot}{ordipointlabel}(x, ...)
}

\arguments{
  \item{x}{For \code{ordipointlabel} a result object from an
    ordination function or an ordination plot (possibly in a pipe). For
    \code{plot.ordipointlabel} an object from \code{ordipointlabel}.}
  \item{display}{Scores displayed in the plot. The default is to show
    \code{"sites"} and \code{"species"} that are available to many
    ordination methods, but there can be only one set or more than two
    set of scores to \code{display}. }
  \item{choices}{Axes shown. }
  \item{col, pch, font, cex}{Colours, point types, font style and
     character expansion for each kind of scores displayed in the
     plot. These should be vectors of the same length as the number of
     items in \code{display}, or if there is only one \code{display}
     they can be a vector of the length of number items.}
  \item{add}{ Add to an existing plot. Default is \code{add = TRUE} when
     the function is used in a pipe, and \code{FALSE} usually. }

   \item{labels}{Labels used in graph. Species (variable) and SU (row)
    names are used if this is missing. Labels must be given in one
    vector for all scores of \code{display}. Function \code{labels} can
    extract the current name from a saved \code{ordipointlabel}
    object. If \code{select} is used, labels are given only for the
    selected items.}
  
  \item{bg}{Background colour for labels. If this is given, texts is
    drawn over non-transparent background. Either a single colour or
    vector of colours for each \code{display}, or with one display, for
    each label.}
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
  \code{method = "SANN"}) to optimize the locations of the text labels
  to the points. There are eight possible locations: up, down, two sides
  and four corners. There is a weak preference to text away from zero,
  and a weak avoidance of corners. The locations and goodness of
  solution varies between runs, and there is no guarantee of finding the
  global optimum, or the same text locations twice. The optimization can
  take a long time in difficult cases with a high number of potential
  overlaps. Several sets of scores can be displayed in one plot.

  The function can be used in a pipe where the first command is an
  ordination \code{plot} command with \code{type = "n"} or to add
  points and lablels to save \pkg{vegan} ordination plot object. See
  examples.

}
\value{
  The function returns invisibly an object of class
  \code{ordipointlabel} with items \code{xy} for coordinates of
  points, \code{labels} for coordinates of labels, items \code{pch},
  \code{cex} and \code{font} for graphical parameters of each point or
  label. In addition, it returns the result of \code{\link{optim}} as
  an attribute \code{"optim"}. The unit of overlap is the area
  of character \code{"m"}, and with varying graphical parameters the
  smallest alternative.

  There is a \code{plot} method based on \code{orditkplot} but it
  does not alter or reset the graphical parameters via \code{par}.
  
  The result object from \code{ordipointlabel} is similar as the
  \code{orditkplot} object of \CRANpkg{vegan3d} package, and it may
  be possible to further edit the result object with
  \code{orditkplot}, but for good results it is
  necessary that the points span the whole horizontal axis without empty
  margins.
}

\seealso{The function is invoked for one set of scores (one
  \code{display}) from \code{\link{text.ordiplot}} and
  \code{\link{plot.cca}} with argument \code{optimize = TRUE}.}

\author{ Jari Oksanen }

\examples{
data(dune, dune.env)
ord <- cca(dune)
ordipointlabel(ord)

## Use in a pipe: optimize species, sites & centroids together
ord <- cca(dune ~ Management + Moisture, dune.env)
plot(ord, scaling = "symmetric", type = "n") |>
  ordipointlabel(c("sites","species","centroids"), cex=c(0.7,0.7,1),
      col = c("black","red","blue"), font = c(1,3,1), pch=c(1,3,4), xpd=TRUE) |>
  text("biplot", col = "blue", bg = "white", cex=1)
}
\keyword{ hplot }
\keyword{ aplot }

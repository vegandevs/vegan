\name{ordiArrowTextXY}
\alias{ordiArrowMul}
\alias{ordiArrowTextXY}

\title{Support Functions for Drawing Vectors}
\description{
  Support functions to assist with drawing of vectors (arrows) on
  ordination plots. \code{ordiArrowMul} finds the multiplier for the
  coordinates of the head of the vector such that they accupy
  \code{fill} proportion of the plot region. \code{ordiArrowTextXY}
  finds coordinates for the locations of \code{labels} to be drawn just
  beyond the head of the vector.
}
\usage{
ordiArrowTextXY(x, labels, display, choices = c(1,2),
                rescale = TRUE, fill = 0.75, \ldots)
ordiArrowMul(x, at = c(0,0), fill = 0.75,
             display, choices = c(1,2), \ldots)
}

\arguments{
  \item{x}{An R object, from which \code{\link{scores}} can determine
    suitable ordination scores or an object created by
    \code{\link{envfit}}, or a two-column matrix of coordinates of arrow
    heads on the two plot axes.}
    
  \item{labels}{Change plotting labels. A character vector of labels for
    which label coordinates are sought. If not supplied, these will be
    determined from the row names of \code{x}, or \code{scores(x, ...)}
    if required. If either of these are not defined, suitable labels
    will be generated.}
  
  \item{display}{a character string known to \code{\link{scores}} or one
    of its methods which indicates the type of scores to extract. In
    fitting functions these are ordinary site scores or linear
    combination scores (\code{"lc"}) in constrained ordination
    (\code{\link{cca}}, \code{\link{rda}}, \code{\link{capscale}}). If
    \code{x} was created by \code{envfit} then \code{display} can not be
    set by the user and takes the value \code{"vectors"}. Ignored if
    \code{x} is a matrix.}
    
  \item{choices}{Axes to be plotted.}

  \item{rescale}{logical; should the coordinates in or extracted from
    \code{x} be rescaled to fill \code{fill} proportion of the plot
    region? The default is to always rescale the coordinates as this is
    usually desired for objects \code{x} from which coordinates are
    retrieved. If supplying \code{x} a 2-column matrix that has already
    been rescaled, then set this to \code{FALSE}.}
    
  \item{fill}{numeric; the proportion of the plot to fill by the span of
    the arrows.}
  
  \item{at}{The origin of fitted arrows in the plot.  If you plot arrows
    in other places than origin, you probably have to specify
    \code{arrrow.mul}.}
  
  \item{...}{Parameters passed to \code{\link{scores}}, and
    \code{\link{strwidth}} and \code{\link{strheight}}.}
}

\details{
  \code{ordiArrowMul} finds a multiplier to scale a bunch of
  arrows to fill an ordination plot, and \code{ordiArrowTextXY} finds
  the coordinates for labels of these arrows. NB.,
  \code{ordiArrowTextXY} does not draw labels; it simply returns
  coordinates at which the labels should be drawn for use with another
  function, such as \code{\link{text}}.
}

\value{
  For \code{ordiArrowTextXY}, a 2-column matrix of coordinates for the
  label centres in the coordinate system of the currently active
  plotting device.

  For \code{ordiArrowMul}, a length-1 vector containing the scaling
  factor.
}

\author{Jari Oksanen, with modifications by Gavin L. Simpson}

\examples{
    ## Scale arrows by hand to fill 80% of the plot
    ## Biplot arrows by hand
    data(varespec, varechem)
    ord <- cca(varespec ~ Al + P + K, varechem)
    plot(ord, display = c("species","sites"))

    ## biplot scores
    bip <- scores(ord, choices = 1:2, display = "bp")

    ## scaling factor for arrows to fill 80% of plot
    (mul <- ordiArrowMul(bip, fill = 0.8))
    bip.scl <- bip * mul                    # Scale the biplot scores
    labs <- rownames(bip)                   # Arrow labels

    ## calculate coordinate of labels for arrows
    (bip.lab <- ordiArrowTextXY(bip.scl, rescale = FALSE, labels = labs))

    ## draw arrows and text labels
    arrows(0, 0, bip.scl[,1], bip.scl[,2], length = 0.1)
    text(bip.lab, labels = labs)

    ## Handling of ordination objects directly
    mul2 <- ordiArrowMul(ord, display = "bp", fill = 0.8)
    stopifnot(all.equal(mul, mul2))
}

\keyword{utilities}
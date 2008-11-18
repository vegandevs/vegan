\name{linestack}
\alias{linestack}

\title{Plots One-dimensional Diagrams without Overwriting Labels }
\description{
  Function \code{linestack} plots vertical one-dimensional plots for
  numeric vectors.  The plots are always labelled, but the labels
  are moved vertically to avoid overwriting.
}
\usage{
linestack(x, labels, cex = 0.8, side = "right", hoff = 2, air = 1.1,
    at = 0, add = FALSE, axis = FALSE, ...)
}

\arguments{
  \item{x}{Numeric vector to be plotted. }
  \item{labels}{Text labels used instead of default (names of \code{x}).}
  \item{cex}{Size of the labels. }
  \item{side}{Put labels to the \code{"right"} or
    \code{"left"} of the axis. }
  \item{hoff}{Distance from the vertical axis to the label in units of
    the width of letter \dQuote{m}. }
  \item{air}{Multiplier to string height to leave empty space between labels. }
  \item{at}{Position of plot in horizontal axis. }
  \item{add}{Add to an existing plot. }
  \item{axis}{Add axis to the plot. }
  \item{\dots}{Other graphical parameters to labels.}
}
\value{
  The function returns invisibly the shifted positions of labels in
  user coordinates.
}
\author{Jari Oksanen }
\note{ The function always draws labelled diagrams.  If you want to have
  unlabelled diagrams, you can use, e.g., \code{\link{plot}},
  \code{\link{stripchart}} or \code{\link{rug}}.
}

\examples{
## First DCA axis
data(dune)
ord <- decorana(dune)
linestack(scores(ord, choices=1, display="sp"))
linestack(scores(ord, choices=1, display="si"), side="left", add=TRUE)
title(main="DCA axis 1")
}
\keyword{ hplot }
\keyword{ aplot }

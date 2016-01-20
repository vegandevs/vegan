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
  \item{labels}{Labels used instead of default (names of \code{x}). May
    be expressions to be drawn with \code{\link{plotmath}}.}
  \item{cex}{Size of the labels. }
  \item{side}{Put labels to the \code{"right"} or \code{"left"} of the
    axis.}
  \item{hoff}{Distance from the vertical axis to the label in units of
    the width of letter \dQuote{m}. }
  \item{air}{Multiplier to string height to leave empty space between
    labels.}
  \item{at}{Position of plot in horizontal axis. }
  \item{add}{Add to an existing plot. }
  \item{axis}{Add axis to the plot. }
  \item{\dots}{Other graphical parameters to labels.}
}
\value{
  The function returns invisibly the shifted positions of labels in
  user coordinates.
}
\author{Jari Oksanen with modifications by Gavin L. Simpson}
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

## Expressions as labels
N <- 10					# Number of sites
df <- data.frame(Ca = rlnorm(N, 2), NO3 = rlnorm(N, 4),
                 SO4 = rlnorm(N, 10), K = rlnorm(N, 3))
ord <- rda(df, scale = TRUE)
### vector of expressions for labels
labs <- expression(Ca^{2+phantom()},
                   NO[3]^{-phantom()},
                   SO[4]^{-2},
                   K^{+phantom()})
scl <- "sites"
linestack(scores(ord, choices = 1, display = "species", scaling = scl),
          labels = labs, air = 2)
linestack(scores(ord, choices = 1, display = "site", scaling = scl),
          side = "left", add = TRUE)
title(main = "PCA axis 1")
}
\keyword{ hplot }
\keyword{ aplot }

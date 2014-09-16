\name{orditkplot}
\alias{orditkplot}
\alias{plot.orditkplot}
\alias{scores.orditkplot}
\alias{points.orditkplot}
\alias{text.orditkplot}

\title{ Ordination Plot with Movable Labels }
\description{
  Function \code{orditkplot} produces an editable ordination plot with
  points and labels. The labels can be moved with mouse, and the edited
  plot can be saved as an encapsulated postscript file or exported via \R 
  \code{plot} function to other graphical formats, or saved in the \R session 
  for further processing.
}
\usage{
orditkplot(x, display = "species", choices = 1:2, width, xlim, ylim, 
   tcex = 0.8, tcol, pch = 1,  pcol, pbg, pcex = 0.7, labels,  ...)
\method{plot}{orditkplot}(x, ...)
\method{points}{orditkplot}(x, ...)
\method{text}{orditkplot}(x, ...)
\method{scores}{orditkplot}(x, display, ...)
}

\arguments{ 
  \item{x}{An ordination result or any other object that
    \code{\link{scores}} can handle, or for the \code{plot} function the
    object dumped from the interactive \code{orditkplot} session.}
  \item{display}{Type of \code{\link{scores}} displayed. For ordination
    scores this typically is either \code{"species"} or \code{"sites"},
    and for \code{orditkplot} result it is either \code{"points"} or
    \code{"labels"}.}
  \item{choices}{Axes displayed.}
  \item{width}{Width of the plot in inches; defaults to the current
    width of the graphical device. }
  \item{xlim, ylim}{x and y limits for plots: points outside these
    limits will be completely removed.}
  \item{tcex}{Character expansion for text labels.}
  \item{tcol}{Colour of text labels.}
  \item{pch, pcol, pbg}{Point type and outline and fill colours. 
    Defaults \code{pcol="black"}  and \code{pbg="transparent"}. 
   Argument \code{pbg} has an effect
    only in filled plotting characters \code{pch = 21} to \code{25}.  } 
  \item{pcex}{Expansion factor for point size.}  
  \item{labels}{Labels used instead of row names.}
  \item{\dots}{Other arguments passed to the function. These can be
    graphical parameters (see \code{\link{par}}) used in the plot, or
    extra arguments to \code{\link{scores}}. These arguments are
    ignored in \code{plot}, but honoured in \code{text} and \code{points}. }  
}

\details{ Function \code{orditkplot} uses \pkg{tcltk} package to draw
  Tcl/Tk based ordination graphics with points and labels. The function
  opens an editable canvas with fixed points, but the labels can be
  dragged with mouse to better positions or edited. In addition, it is
  possible to zoom to a part of the graph.

  The function knows the following mouse operations:
  \itemize{
    \item \strong{Left mouse button} can be used to move labels to better
    positions. A line will connect a label to the corresponding point.
    \item \strong{Double clicking left mouse button} opens a window where the
    label can be edited. After editing the label, hit the Return key.
    \item \strong{Right mouse button} (or alternatively, Shift-Mouse button with
    one-button mouse) can be used for zooming to a part of the
    graph. Keeping the mouse button down and dragging will draw a box
    of the zoomed area, and after releasing the button, a new plot window
    will be created (this is still preliminary: all arguments are not
    passed to the new plot).
  }

  In addition there are buttons
  for the following tasks: \strong{Copy to EPS} copies the current plot
  to an encapsulated postscript (eps) file using standard Tcl/Tk
  utilities. The faithfulness of this copy is system dependent. Button
  \strong{Export plot} uses \code{plot.orditkplot} function to redraw
  the plot into graphical file formats. Depending on the system, the
  following graphical formats may be available: eps, pdf, png, jpeg or
  bmp. The file type is deduced from the file suffix or the selection of
  the file type in the dialogue box. Alternatively, the same dialogue can be
  used to save the plot to an editable \code{\link{xfig}} file. Button
  \strong{Dump to R} writes the edited coordinates of labels and points
  to the \R session for further processing, and the
  \code{plot.orditkplot} function can be used to display the
  results. For faithful replication of the plot, the graph must have
  similar dimensions as the \code{orditkplot} canvas had originally. The
  \code{plot} function cannot be configured, but it uses the same
  settings as the original Tcl/Tk plot. However, \code{points} and
  \code{text} functions are fully configurable, and unaware of the
  original Tcl/Tk plot settings (probably you must set \code{cex} at
  least to get a decent plot). Finally, button \strong{Dismiss} closes
  the window.

  The produced plot will have equal aspect ratio. The width of the
  horizontal axis is fixed, but vertical axes will be scaled to needed
  height, and you can use scrollbar to move vertically if the whole
  canvas does not fit the window. If you use dumped labels in ordinary
  \R plots, your plot must have the same dimensions as the
  \code{orditkplot} canvas to have identical location of the labels.

  The function only displays one set of scores. However, you can use
  \code{\link{ordipointlabel}} to produce a result object that has
  different points and text types for several sets of scores and this
  can further edited with \code{orditkplot}. For a good starting
  solution you need to scale the \code{\link{ordipointlabel}} result so
  that the points span over the whole horizontal axis.

  The plot is a Tcl/Tk canvas, but the function tries to replicate
  standard graphical device of the platform, and it honours several
  graphical parameters (see \code{\link{par}}).  Many of the graphical
  parameters can be given on the command line, and they will be passed
  to the function without influencing other graphical devices in \R.
  At the moment, the
  following graphical parameters are honoured: 
  \code{pch} \code{bg},  \code{cex},
  \code{cex.axis}, \code{cex.lab}, \code{col} (for labels),
  \code{col.axis}, \code{col.lab}, \code{family} (for font faces),
  \code{fg}, \code{font}, \code{font.axis}, \code{font.lab},
  \code{lheight}, \code{lwd} (for the box), \code{mar}, \code{mex},
  \code{mgp}, \code{ps}, \code{tcl}. These can be set with
  \code{\link{par}}, and they also will influence other plots similarly.
  
  The \code{\link[tcltk]{tkcanvas}} text cannot be rotated, and
  therefore vertical axis is not labelled, and \code{las}
  \code{\link{par}}ameter will not be honoured in the Tcl/Tk plot, but
  it will be honoured in the exported \R plots and in
  \code{plot.orditkplot}.  }

\value{ Function returns nothing useful directly, but you can save the
  edited graph to a file or dump the edited positions to an \R session
  for further processing and plotting.  
}

\author{ Jari Oksanen }
\note{
  You need \pkg{tcltk} package and \R must have been configured with
  \code{\link{capabilities}} for \code{tcltk} when building the binary.
  Depending on your OS, you may need to start X11 and set the display
  before loading \pkg{tcltk} and starting the function (for instance,
  with \code{Sys.setenv("DISPLAY"=":0")}). See
  \code{\link[tcltk]{tcltk-package}}. 

}
\seealso{ Function \code{\link{ordipointlabel}} is an automatic
  procedure with similar goals of avoiding overplotting. 	   
  See \code{\link{ordiplot}}, \code{\link{plot.cca}},
  \code{\link[vegan3d]{ordirgl}} and \code{\link{orditorp}} for alternative
  ordination plots, and \code{\link{scores}} for extracting ordination
  scores.  }
\examples{
## The example needs user interaction and is not executed directly.
## It should work when pasted to the window.
\dontrun{
data(varespec)
ord <- cca(varespec)
## Do something with the graph and end by clicking "Dismiss"
orditkplot(ord, mar = c(4,4,1,1)+.1, font=3)
## Use ordipointlabel to produce a plot that has both species and site
## scores in different colors and plotting symbols
pl <- ordipointlabel(ord)
orditkplot(pl)
}
}
\keyword{ iplot }
\keyword{ dynamic }

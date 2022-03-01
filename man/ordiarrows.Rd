\name{ordiarrows}
\alias{ordiarrows}
\alias{ordisegments}
\alias{ordigrid}

\title{Add Arrows and Line Segments to Ordination Diagrams}

\description{ Functions to add arrows, line segments, regular grids of
  points. The ordination diagrams can be produced by \code{vegan}
  \code{\link{plot.cca}}, \code{\link{plot.decorana}} or
  \code{\link{ordiplot}}.  }

\usage{
ordiarrows(ord, groups, levels, replicates, order.by, display = "sites",
         col = 1, show.groups, startmark, label = FALSE, length = 0.1, ...)
ordisegments(ord, groups, levels, replicates, order.by, display = "sites",
         col = 1, show.groups, label = FALSE, ...)
ordigrid(ord, levels, replicates, display = "sites",  lty = c(1,1), 
         col = c(1,1), lwd = c(1,1), ...)
}

\arguments{
  \item{ord}{An ordination object or an \code{\link{ordiplot}} object. }
  \item{groups}{Factor giving the groups for which the graphical item is
    drawn. }
  \item{levels, replicates}{Alternatively, regular
    groups can be defined with arguments \code{levels} and
    \code{replicates}, where \code{levels} gives the number of groups,
    and \code{replicates} the number of successive items at the same
    group.}
  \item{order.by}{Order points by increasing order of this variable
    within \code{groups}. Reverse sign of the variable for decreasing
    ordering.} 
  \item{display}{Item to displayed. }
  \item{show.groups}{Show only given groups. This can be a vector, or
    \code{TRUE} if you want to show items for which condition is
    \code{TRUE}. This argument makes it possible to use different
    colours and line types for groups. The default is to show all groups. }
  \item{label}{Label the \code{groups} by their names. In
    \code{ordiellipse}, \code{ordihull} and \code{ordispider} the the
    group name is in the centroid of the object, in \code{ordiarrows}
    in the start of the arrow, and in \code{ordisegments} at both
    ends. \code{ordiellipse} and \code{ordihull} use standard
    \code{\link{text}}, and others use \code{\link{ordilabel}}.}
  
  \item{startmark}{plotting character used to mark the first item. The
    default is to use no mark, and for instance, \code{startmark = 1}
    will draw a circle.  For other plotting characters, see \code{pch}
    in \code{\link{points}}. }

  \item{col}{Colour of lines, \code{label} borders and
    \code{startmark} in \code{ordiarrows} and
    \code{ordisegments}. This can be a vector recycled for
    \code{groups}. In \code{ordigrid} it can be a vector of length 2
    used for \code{levels} and \code{replicates}.}

  \item{length}{Length of edges of the arrow head (in inches).}

  \item{lty, lwd}{Line type, line width used for 
    \code{level}s and \code{replicate}s in \code{ordigrid}.}

  \item{\dots}{Parameters passed to graphical functions such as
    \code{\link{lines}}, \code{\link{segments}}, \code{\link{arrows}},
    or to \code{\link{scores}} to select axes and scaling etc. }
}

\details{ Function \code{ordiarrows} draws \code{\link{arrows}} and
  \code{ordisegments} draws line \code{\link{segments}} between
  successive items in the groups. Function \code{ordigrid} draws line
  \code{\link{segments}} both within the groups and for the
  corresponding items among the groups.

}

\note{These functions add graphical items to ordination graph: You must
  draw a graph first.
}
\author{ Jari Oksanen }

\seealso{The functions pass parameters to basic graphical functions, and
  you may wish to change the default values in \code{\link{arrows}},
  \code{\link{lines}} and \code{\link{segments}}. You can pass
  parameters to \code{\link{scores}} as well. 
}

\examples{
example(pyrifos)
mod <- rda(pyrifos)
plot(mod, type = "n")
## Annual succession by ditches, colour by dose
ordiarrows(mod, ditch, label = TRUE, col = as.numeric(dose))
legend("topright", levels(dose), lty=1, col=1:5, title="Dose")
## Show only control and highest Pyrifos treatment
plot(mod, type = "n")
ordiarrows(mod, ditch, label = TRUE, 
   show.groups = c("2", "3", "5", "11"))
ordiarrows(mod, ditch, label = TRUE, show = c("6", "9"),
   col = 2)
legend("topright", c("Control", "Pyrifos 44"), lty = 1, col = c(1,2))
}
\keyword{aplot }

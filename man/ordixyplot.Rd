\name{ordixyplot}
\alias{ordixyplot}
\alias{ordisplom}
\alias{ordicloud}
\alias{panel.ordi}
\alias{panel.ordiarrows}
\alias{panel.ordi3d}
\alias{prepanel.ordi3d}
\alias{ordilattice.getEnvfit}

\title{ Trellis (Lattice) Plots for Ordination }
\description{
  Functions \code{ordicloud}, \code{ordisplom} and \code{ordixyplot}
  provide an interface to plot ordination results using Trellis
  functions \code{\link[lattice]{cloud}}, \code{\link[lattice]{splom}}
  and \code{\link[lattice]{xyplot}} in package \pkg{lattice}. 
}
\usage{
ordixyplot(x, data = NULL, formula, display = "sites", choices = 1:3,
    panel = "panel.ordi", aspect = "iso", envfit,
    type = c("p", "biplot"), ...)
ordisplom(x, data=NULL, formula = NULL,  display = "sites", choices = 1:3,
    panel = "panel.ordi", type = "p",  ...)
ordicloud(x, data = NULL, formula, display = "sites", choices = 1:3, 
    panel = "panel.ordi3d", prepanel = "prepanel.ordi3d",  ...)
}
\arguments{
  \item{x}{ An ordination result that \code{\link{scores}} knows: any
    ordination result in \pkg{vegan} and many others. }
  \item{data}{ Optional data to amend ordination results. The ordination
    results are found from \code{x}, but you may give here data for other
    variables needed in plots. Typically these are environmental data.}
  \item{formula}{ Formula to define the plots. A default formula will be
    used if this is omitted. The
    ordination axes must be called by the same names as in the
    ordination results (and these names vary among methods). In
    \code{ordisplom}, special character \code{.} refers to the
    ordination result. }
  \item{display}{ The kind of scores: an argument passed to
    \code{\link{scores}}. }
  \item{choices}{ The axes selected: an argument passed to
    \code{\link{scores}}. } 
  \item{panel, prepanel}{ The names of the panel and prepanel
    functions. } 
  \item{aspect}{The aspect of the plot (passed to the \pkg{lattice}
    function).}
  \item{envfit}{Result of \code{\link{envfit}} function displayed in
    \code{ordixyplot}. Please note that this needs same \code{choices}
    as \code{ordixyplot}.}
  \item{type}{The type of plot. This knows the same alternatives as
    \code{\link[lattice]{panel.xyplot}}. In addition \code{ordixyplot}
    has alternatives \code{"biplot"} and \code{"arrows"}. The first
    displays fitted vectors and 
    factor centroids of \code{envfit}, or in constrained ordination, the
    biplot arrows and factor centroids if \code{envfit} is not
    given. The second (\code{type = "arrows"}) is a trellis variant of
    \code{\link{ordiarrows}} and draws arrows by \code{groups}. The
    line parameters are controlled by
    \code{\link[lattice]{trellis.par.set}} for \code{superpose.line},
    and the user can set \code{length}, \code{angle} and \code{ends}
    parameters of \code{\link[lattice]{panel.arrows}}.}
  \item{\dots}{ Arguments passed to \code{\link{scores}} methods or
    \pkg{lattice} functions.}
}
\details{
  The functions provide an interface to the corresponding \pkg{lattice}
  functions.  All graphical parameters are passed to the \pkg{lattice}
  function so that these graphs are extremely configurable. See
  \code{\link[lattice]{Lattice}} and \code{\link[lattice]{xyplot}},
    \code{\link[lattice]{splom}} and \code{\link[lattice]{cloud}} for
  details, usage and possibilities.

  The argument \code{x} must always be an ordination result. The scores
  are extracted with \pkg{vegan} function \code{\link{scores}} so that
  these functions work with all \pkg{vegan} ordinations and many others.

  The \code{formula} is used to define the models. All functions have
  simple default formulae which are used if \code{formula} is missing. 
  If formula is omitted in \code{ordisplom} it
  produces a pairs plot of ordination axes and variables in
  \code{data}. If \code{formula} is given, ordination results must be
  referred to as \code{.} and other variables by their names. In other
  functions, the formula must use the names of ordination scores and names
  of \code{data}.

  The ordination scores are found from \code{x}, and \code{data} is
  optional. The \code{data} should contain other variables than
  ordination scores to be used in plots. Typically, they are
  environmental variables (typically factors) to define panels or plot
  symbols.

  The proper work is done by the panel function. The layout can be
  changed by defining own panel functions. See
  \code{\link[lattice]{panel.xyplot}},
  \code{\link[lattice]{panel.splom}} and
  \code{\link[lattice]{panel.cloud}} for details and survey of
  possibilities.

  Ordination graphics should always be isometric: same scale should be
  used in all axes. This is controlled (and can be changed) with
  argument \code{aspect} in \code{ordixyplot}. In \code{ordicloud} the
  isometric scaling is defined in \code{panel} and \code{prepanel}
  functions. You must replace these functions if you want to have
  non-isometric scaling of graphs. You cannot select isometric scaling
  in \code{ordisplom}.
}

\value{
  The function return \code{\link[lattice]{Lattice}} objects of class
  \code{"trellis"}.   
}
\author{Jari Oksanen }
 
\seealso{
  \code{\link[lattice]{Lattice}},
  \code{\link[lattice]{xyplot}},
  \code{\link[lattice]{splom}},
  \code{\link[lattice]{cloud}},
  \code{\link[lattice]{panel.splom}},
  \code{\link[lattice]{panel.cloud}}
}
\examples{
data(dune)
data(dune.env)
ord <- cca(dune)
## Pairs plots
ordisplom(ord)
ordisplom(ord, data=dune.env, choices=1:2)
ordisplom(ord, data=dune.env, form = ~ . | Management, groups=Manure)
## Scatter plot
ordixyplot(ord, data=dune.env, form = CA1 ~ CA2 | Management,
  groups=Manure)
## Choose a different scaling
ordixyplot(ord, scaling = 3)
## ... Slices of third axis
ordixyplot(ord, form = CA1 ~ CA2 | equal.count(CA3, 4), type = c("g","p"))
## Display environemntal variables
ordixyplot(ord, envfit = envfit(ord ~ Management + A1, dune.env, choices=1:3))
## 3D Scatter plots
ordicloud(ord, form = CA2 ~ CA3*CA1, groups = Manure, data = dune.env)
ordicloud(ord, form = CA2 ~ CA3*CA1 | Management, groups = Manure,
   data = dune.env, auto.key = TRUE, type = c("p","h"))
}

\keyword{ hplot }


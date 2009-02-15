\name{permCheck}
\alias{permCheck}
\alias{numPerms}
\alias{print.permCheck}
\alias{print.summary.permCheck}
\alias{summary.permCheck}
\alias{getNumObs}
\alias{getNumObs.default}
\alias{getNumObs.integer}
\alias{getNumObs.numeric}
\alias{allPerms}
\alias{print.allPerms}
\alias{summary.allPerms}
\alias{print.summary.allPerms}
\alias{permuplot}

\title{Utility functions for permutation schemes}
\description{
  \code{permCheck} provides checking of permutation schemes for
  validity. \code{numPerms} calculates the maximum number of
  permutations possible under the current permutation
  scheme. \code{allPerms} enumerates all possible permutations for the
  given scheme. \code{getNumObs} is a utility function to return the
  number of observations for a range of R and ordination
  objects. \code{permuplot} produces a graphical representation of the
  selected permutation design.
}
\usage{
permCheck(object, control = permControl(), make.all = TRUE)

\method{summary}{permCheck}(object, \dots)

numPerms(object, control = permControl())

allPerms(n, control = permControl(), max = 9999,
         observed = FALSE)

\method{summary}{allPerms}(object, \dots)

getNumObs(object, \dots)

\method{getNumObs}{default}(object, \dots)

\method{getNumObs}{numeric}(object, \dots)

\method{getNumObs}{integer}(object, \dots)

permuplot(n, control = permControl(), col = par("col"),
          hcol = "red", shade = "lightgrey", xlim = NULL, ylim = NULL,
          inset = 0.1, main = NULL, sub = NULL, ann = par("ann"),
          cex = par("cex"), \dots)
}

\arguments{
  \item{object}{an R object. Specifically, for \code{getNumObs} any
    object handled by \code{\link{scores}}, data frames, matrices, and
    numeric and integer  vectors. See Details for a complete
    description, especially for \code{numPerms}. For
    \code{\link{summary.permCheck}} an object of class
    \code{"permCheck"}. For \code{\link{summary.allPerms}} an object of
    class \code{"allPerms"}.} 
  \item{control}{a list of control values describing properties of the
    permutation design, as returned by a call to
    \code{\link{permControl}}.}
  \item{make.all}{logical; should \code{permCheck} generate all
    possible permutations? Useful if want to check permutation design
    but not produce the matrix of all permutations.}
  \item{n}{the number of observations or an 'object' from which the
    number of observations can be determined via \code{getNumObs}.}
  \item{max}{the maximum number of permutations, below which complete
    enumeration will be attempted. See Details.}
  \item{observed}{logical, should the observed ordering of samples be
    returned as part of the complete enumeration? Default is
    \code{FALSE} to facilitate usage in higher level functions.}
  \item{col, xlim, ylim, main, sub, ann, cex}{Graphical parameters.}
  \item{hcol}{Colour to use for highlighting observations and the border
    colour of the polygons drawn when \code{type = "strata"}.}
  \item{shade}{The polygon shading colour (passed to argument \code{col}
    of function \code{\link{polygon}}) when \code{type = "strata"}.}
  \item{inset}{Proportion of range of x and y coordinates to add to the
    plot x and y limits. Used to create a bit of extra space around the
    margin of each plot.}
  \item{\dots}{arguments to other methods. For \code{permuplot}
    graphical parameters can be passed to plotting functions, though
    note that not all parameters will be accepted gracefully at the
    moment.}
}
\details{
  \code{permCheck}, \code{allPerms}, \code{numPerms} and
  \code{permuplot} are utility functions for working with the new
  permutation schemes available in \code{\link{permuted.index2}}.

  \code{permCheck} is used to check the current permutation schemes
  against the object to which it will be applied. It calculates the
  maximum number of possible permutations for the number of observations
  in \code{object} and the permutation scheme described by
  \code{control}. The returned object contains component \code{control},
  an object of class \code{"permControl"} suitably modified if
  \code{permCheck} identifies a problem.

  The main problem is requesting more permutations than possible with
  the number of observations and the permutation design. In such cases,
  \code{nperm} is reduced to equal the number of possible permutations,
  and complete enumeration of all permutations is turned on
  (\code{control\$complete} is set to \code{TRUE}). 

  Alternatively, if the number of possible permutations is low, and less
  than \code{control\$minperm}, it is better to enumerate all possible
  permutations, and as such complete enumeration of all permutations is
  turned  on (\code{control\$complete} is set to \code{TRUE}).

  Function \code{numPerms} returns the number of permutations for the
  passed \code{object} and the selected permutation
  scheme. \code{object} can be one of a data frame, matrix, an object
  for which a scores method exists, or a numeric or integer vector. In
  the case of a numeric or integer vector, a vector of length 1 can be
  used and it will be expanded to a vector of length \code{object}
  (i.e., \code{1:object}) before computing the number of
  permutations. As such, \code{object} can be the number of observations
  not just the object containing the observations.

  Function \code{allPerms} enumerates all possible permutations for the
  number of observations and the selected permutation scheme. It has
  \code{\link{print}} and \code{\link{summary}} methods. \code{allPerms}
  returns a matrix containing all possible permutations, possibly
  containing the observed ordering (if argument \code{observed} is
  \code{TRUE}). The rows of this matrix are the various permutations and
  the columns reflect the number of samples.

  With free permutation designs, and restricted permutation schemes with
  large numbers of observations, there are a potentially huge number of
  possible permutations of the samples. It would be inefficient, not to
  mention incredibly time consuming, to enumerate them all. Storing all
  possible permutations would also become problematic in such cases. To
  control this and guard against trying to evaluate too large a number
  of permutations, if the number of possible permutations is larger than
  \code{max}, \code{allPerms} exits with an error.

  Function \code{getNumObs} is a simple generic function to return the
  number of observations in a range of R objects. The default method
  will work for any object for which a \code{\link{scores}} method
  exists. This includes matrices and data frames, as well as specific
  methods for numeric or integer vectors.

  \code{permuplot} is a graphical utility function, which produces a
  graphical representation of a permutation design. It takes the number
  of observations and an object returned by \code{\link{permControl}} as
  arguments and produces a plot on the currently active device. If
  strata are present in the design, the plotting region is split into
  sufficient plotting regions (one for each stratum), and the design in
  each stratum plotted.

  Free permutation designs are represented by plotting the observation
  number at random x and y coordinates. Series designs (time series or
  line transects) are represented by plotting the observation numbers
  comprising the series in a circle and the start of the permuted series
  is highlighted using colour \code{hcol}. Grid designs are drawn on a
  regular grid and the top left observation in the original grid is
  highlighted using colour \code{hcol}. Note the ordering used is R's
  standard ordering for matrices - columns are filled first.
}
\value{
  For \code{permCheck} a list containing the maximum number of
  permutations possible and an object of class
  \code{"\link{permControl}"}.

  For \code{allPerms}, and object of class \code{"allPerms"}, a matrix
  whose rows are the set of all possible permutations for the supplies
  number of observations and permutation scheme selected. The matrix has
  two additional attributes \code{control} and
  \code{observed}. Attribute \code{control} contains the argument
  \code{control} (possibly updated via \code{permCheck}). Attribute 
  \code{observed} contains argument \code{observed}.

  For \code{numPerms}, the (numeric) number of possible permutations.

  For \code{getNumObs}, the (numeric) number of observations in
  \code{object}.

  For \code{permuplot}, a plot on the currently active device.
}
%\references{
%}
\note{
  In general, mirroring \code{"series"} or \code{"grid"} designs doubles
  or quadruples, respectively,the number of permutations without
  mirroring (within levels of strata if present). This is \strong{not}
  true in two special cases:
  \enumerate{
    \item In \code{"grid"} designs where the number of columns is equal
    to 2, and
    \item In \code{"series"} designs where the number of observations in
    a series is equal to 2.
  }

  For example, with 2 observations there are 2 permutations for
  \code{"series"} designs:
  \enumerate{
    \item 1-2, and
    \item 2-1.
  }
  If these two permutations were mirrored, we would have:
  \enumerate{
    \item 2-1, and
    \item 1-2.
  }
  It is immediately clear that this is the same set of permutations
  without mirroring (if one reorders the rows). A similar situation
  arises in \code{"grid"} designs where the number of \strong{columns}
  per \emph{grid} is equal to 2. Note that the number of rows per
  \emph{grid} is not an issue here.
}
\author{Gavin Simpson}
\seealso{\code{\link{permuted.index2}} and \code{\link{permControl}}.}

\examples{
## use example data from ?pyrifos
example(pyrifos)

## Demonstrate the maximum number of permutations for the pyrifos data
## under a series of permutation schemes

## no restrictions - lots of perms
(check1 <- permCheck(pyrifos, control = permControl(type = "free")))
summary(check1)

## no strata but data are series with no mirroring, so 132 permutations
permCheck(pyrifos, control = permControl(type = "series",
                   mirror = FALSE))

## no strata but data are series with mirroring, so 264 permutations
permCheck(pyrifos, control = permControl(type = "series",
                   mirror = TRUE))

## unrestricted within strata
permCheck(pyrifos, control = permControl(strata = ditch,
                   type = "free"))

## time series within strata, no mirroring
permCheck(pyrifos, control = permControl(strata = ditch,
                   type = "series", mirror = FALSE))

## time series within strata, with mirroring
permCheck(pyrifos, control = permControl(strata = ditch,
                   type = "series", mirror = TRUE))

## time series within strata, no mirroring, same permutation within strata
permCheck(pyrifos, control = permControl(strata = ditch,
                   type = "series", constant = TRUE))

## time series within strata, with mirroring, same permutation within strata
permCheck(pyrifos, control = permControl(strata = ditch,
                   type = "series", mirror = TRUE, constant = TRUE))

## permute strata
permCheck(pyrifos, permControl(strata = ditch, type = "free",
                               permute.strata = TRUE))

## this should also also for arbitrary vectors
vec1 <- permCheck(1:100)
vec2 <- permCheck(1:100, permControl())
all.equal(vec1, vec2)
vec3 <- permCheck(1:100, permControl(type = "series"))
all.equal(100, vec3$n)
vec4 <- permCheck(1:100, permControl(type = "series", mirror = TRUE))
all.equal(vec4$n, 200)

## enumerate all possible permutations
fac <- gl(2,6)
ctrl <- permControl(type = "grid", mirror = FALSE, strata = fac,
                    constant = TRUE, nrow = 3, ncol = 2)
numPerms(1:12, control = ctrl)
(tmp <- allPerms(12, control = ctrl, observed = TRUE))
(tmp2 <- allPerms(12, control = ctrl))
## turn on mirroring
ctrl$mirror <- TRUE
numPerms(1:12, control = ctrl)
(tmp3 <- allPerms(12, control = ctrl, observed = TRUE))
(tmp4 <- allPerms(12, control = ctrl))
## prints out details of the permutation scheme as
## well as the matrix of permutations
summary(tmp)
summary(tmp2)

## different numbers of observations per level of strata
fac <- factor(rep(1:3, times = c(3,2,2)))
## free permutations in levels of strata
numPerms(7, permControl(type = "free", strata = fac))
allPerms(7, permControl(type = "free", strata = fac))
## series permutations in levels of strata
numPerms(7, permControl(type = "series", strata = fac))
allPerms(7, permControl(type = "series", strata = fac))

## allPerms can work with a vector
vec <- c(3,4,5)
allPerms(vec)

## Tests for permuplot
n <- 25
## standard permutation designs
permuplot(n, permControl(type = "free"))
permuplot(n, permControl(type = "series"))
permuplot(n, permControl(type = "grid", nrow = 5, ncol = 5))

## restricted perms with mirroring
permuplot(n, permControl(type = "series", mirror = TRUE))
permuplot(n, permControl(type = "grid", nrow = 5, ncol = 5,
                             mirror = TRUE))

## perms within strata
fac <- gl(6, 20)
control <- permControl(type = "free", strata = fac)
permuplot(120, control = control, cex = 0.8)
control <- permControl(type = "series", strata = fac)
permuplot(120, control = control, cex = 0.8)
fac <- gl(6, 25)
control <- permControl(type = "grid", strata = fac,
                       nrow = 5, ncol = 5)
permuplot(150, control = control, cex = 0.8)

## perms within strata with mirroring
fac <- gl(6, 20)
control <- permControl(type = "series", strata = fac,
                       mirror = TRUE)
permuplot(120, control = control, cex = 0.8)
fac <- gl(6, 25)
control <- permControl(type = "grid", strata = fac,
                       nrow = 5, ncol = 5, mirror = TRUE)
permuplot(150, control = control, cex = 0.8)

## same perms within strata
fac <- gl(6, 20)
control <- permControl(type = "free", strata = fac,
                       constant = TRUE)
permuplot(120, control = control, cex = 0.8)
control <- permControl(type = "series", strata = fac,
                       constant = TRUE)
permuplot(120, control = control, cex = 0.8)
fac <- gl(6, 25)
control <- permControl(type = "grid", strata = fac,
                       nrow = 5, ncol = 5, constant = TRUE)
permuplot(150, control = control, cex = 0.8)

## same perms within strata with mirroring
fac <- gl(6, 20)
control <- permControl(type = "series", strata = fac,
                       mirror = TRUE, constant = TRUE)
permuplot(120, control = control, cex = 0.8)
fac <- gl(6, 25)
control <- permControl(type = "grid", strata = fac,
                       nrow = 5, ncol = 5, mirror = TRUE,
                       constant = TRUE)
permuplot(150, control = control, cex = 0.8)

}
\keyword{ utilities }
\keyword{ design }
\keyword{ methods }
\keyword{ datagen }

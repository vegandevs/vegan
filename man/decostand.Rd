\encoding{UTF-8}
\name{decostand}
\alias{decostand}
\alias{wisconsin}
\alias{decobackstand}

\title{Standardization Methods for Community Ecology}
\description{
The function provides some popular (and effective) standardization
methods for community ecologists.
}
\usage{
decostand(x, method, MARGIN, range.global, logbase = 2, na.rm=FALSE, ...)
wisconsin(x)
decobackstand(x, zap = TRUE)
}

\arguments{
  \item{x}{Community data, a matrix-like object. For
    \code{decobackstand} standardized data. }
  \item{method}{Standardization method. See Details for available options.}
  \item{MARGIN}{Margin, if default is not acceptable. \code{1} = rows,
    and \code{2} = columns of \code{x}.}
  \item{range.global}{Matrix from which the range is found in
    \code{method = "range"}.  This allows using same ranges across
    subsets of data.  The dimensions of \code{MARGIN} must match with
    \code{x}. }
  \item{logbase}{The logarithm base used in \code{method = "log"}.}
  \item{na.rm}{Ignore missing values in row or column standardizations.}
  \item{zap}{Make near-zero values exact zeros to avoid negative
    values and exaggerated estimates of species richness.}
  \item{\dots}{Other arguments to the function (ignored).}
}
\details{
  The function offers following standardization methods for community
  data:
  \itemize{
    \item \code{total}: divide by margin total (default \code{MARGIN = 1}).
    \item \code{max}: divide by margin maximum (default \code{MARGIN = 2}).
    \item \code{frequency}: divide by margin total and multiply by the
    number of non-zero items, so that the average of non-zero entries is
    one (Oksanen 1983; default \code{MARGIN = 2}).
    \item \code{normalize}: make margin sum of squares equal to one (default
    \code{MARGIN = 1}).
    \item \code{range}: standardize values into range 0 \dots 1 (default
    \code{MARGIN = 2}).  If all values are constant, they will be
    transformed to 0.
    
    \item \code{rank, rrank}: \code{rank} replaces abundance values by
    their increasing ranks leaving zeros unchanged, and \code{rrank} is
    similar but uses relative ranks with maximum 1 (default
    \code{MARGIN = 1}). Average ranks are used for tied values.

    \item \code{standardize}: scale \code{x} to zero mean and unit variance
    (default \code{MARGIN = 2}).
    \item \code{pa}: scale \code{x} to presence/absence scale (0/1).
    \item \code{chi.square}: divide by row sums and square root of
    column sums, and adjust for square root of matrix total
    (Legendre & Gallagher 2001). When used with the Euclidean
    distance, the distances should be similar to the
    Chi-square distance used in correspondence analysis. However, the
    results from \code{\link{cmdscale}} would still differ, since
    CA is a weighted ordination method (default \code{MARGIN = 1}).
    \item \code{hellinger}: square root of \code{method = "total"}
    (Legendre & Gallagher 2001).
    \item \code{log}: logarithmic transformation as suggested by
     Anderson et al. (2006): \eqn{\log_b (x) + 1}{log_b (x) + 1} for
     \eqn{x > 0}, where \eqn{b} is the base of the logarithm; zeros are
     left as zeros. Higher bases give less weight to quantities and more
     to presences, and \code{logbase = Inf} gives the presence/absence
     scaling. Please note this is \emph{not} \eqn{\log(x+1)}{log(x+1)}.
     Anderson et al. (2006) suggested this for their (strongly) modified
     Gower distance (implemented as \code{method = "altGower"} in 
     \code{\link{vegdist}}), but the standardization can be used 
     independently of distance indices.

     \item \code{alr}: Additive log ratio ("alr") transformation
     (Aitchison 1986) reduces data skewness and compositionality
     bias. The transformation assumes positive values, pseudocounts can
     be added with the argument \code{pseudocount}. One of the
     rows/columns is a reference that can be given by \code{reference}
     (name of index). The first row/column is used by default
     (\code{reference = 1}).  Note that this transformation drops one
     row or column from the transformed output data. The \code{alr}
     transformation is defined formally as follows:
     \deqn{alr = [log\frac{x_1}{x_D}, ..., log\frac{x_{D-1}}{x_D}]}{%
       alr = [log(x_i/x_D), ..., log(x_{-D}/x_D)]}
     where the denominator sample \eqn{x_D} can be chosen
     arbitrarily. This transformation is often used with pH and other
     chemistry measurenments. It is also commonly used as multinomial
     logistic regression. Default \code{MARGIN = 1} uses row as the
     \code{reference}.

     \item \code{clr}: centered log ratio ("clr") transformation proposed by
     Aitchison (1986) and it is used to reduce data skewness and compositionality bias.
     This transformation has frequent applications in microbial ecology
     (see e.g. Gloor et al., 2017). The clr transformation is defined as:
     \deqn{clr = log\frac{x}{g(x)} = log x - log g(x)}{%
           clr = log(x/g(x)) = log x - log g(x)}     
     where \eqn{x} is a single value, and g(x) is the geometric mean of
     \eqn{x}.
     The method can operate only with positive data;
     a common way to deal with zeroes is to add pseudocount
     (e.g. the smallest positive value in the data), either by
     adding it manually to the input data, or by using the argument
     \code{pseudocount} as in
     \code{decostand(x, method = "clr", pseudocount = 1)}. Adding
     pseudocount will inevitably introduce some bias; see
     the \code{rclr} method for one available solution.

     \item \code{rclr}: robust clr ("rclr") is similar to regular clr
     (see above) but allows data that contains zeroes. This method
     does not use pseudocounts, unlike the standard clr.
     The robust clr (rclr) divides the values by geometric mean
     of the observed features; zero values are kept as zeroes, and not
     taken into account. In high dimensional data,
     the geometric mean of rclr approximates the true
     geometric mean; see e.g. Martino et al. (2019)
     The \code{rclr} transformation is defined formally as follows:
     \deqn{rclr = log\frac{x}{g(x > 0)}}{%
     rclr = log(x/g(x > 0))}
     where \eqn{x} is a single value, and \eqn{g(x > 0)} is the geometric 
     mean of sample-wide values \eqn{x} that are positive (> 0).
  }
  
  Standardization, as contrasted to transformation, means that the
  entries are transformed relative to other entries.

  All methods have a default margin. \code{MARGIN=1} means rows (sites
  in a normal data set) and \code{MARGIN=2} means columns (species in a
  normal data set).

  Command \code{wisconsin} is a shortcut to common Wisconsin double
  standardization where species (\code{MARGIN=2}) are first standardized
  by maxima (\code{max}) and then sites (\code{MARGIN=1}) by
  site totals (\code{tot}).

  Most standardization methods will give nonsense results with
  negative data entries that normally should not occur in the community
  data. If there are empty sites or species (or constant with
  \code{method =  "range"}), many standardization will change these into
  \code{NaN}.

  Function \code{decobackstand} can be used to transform standardized
  data back to original. This is not possible for all standardization
  and may not be implemented to all cases where it would be
  possible. There are round-off errors and back-transformation is not
  exact, and it is wise not to overwrite the original data. With
  \code{zap=TRUE} original zeros should be exact.
}

\value{
  Returns the standardized data frame, and adds an attribute
  \code{"decostand"} giving the name of applied standardization
  \code{"method"} and attribute \code{"parameters"} with appropriate
  transformation parameters.  }

\author{Jari Oksanen, Etienne \enc{Laliberté}{Laliberte}
  (\code{method = "log"}), Leo Lahti (\code{alr}, 
  \code{"clr"} and \code{"rclr"}).}
  
\note{Common transformations can be made with standard \R functions.}

\references{ 

  Aitchison, J. The Statistical Analysis of Compositional Data (1986).
  London, UK: Chapman & Hall.

  Anderson, M.J., Ellingsen, K.E. & McArdle, B.H. (2006) Multivariate
  dispersion as a measure of beta diversity. \emph{Ecology Letters} 
  \strong{9}, 683--693.

  Egozcue, J.J., Pawlowsky-Glahn, V., Mateu-Figueras, G.,
  Barcel'o-Vidal, C. (2003) Isometric logratio transformations for
  compositional data analysis. \emph{Mathematical Geology}
  \strong{35}, 279--300.

  Gloor, G.B., Macklaim, J.M., Pawlowsky-Glahn, V. & Egozcue, J.J. (2017)
  Microbiome Datasets Are Compositional: And This Is Not Optional.
  \emph{Frontiers in Microbiology} \strong{8}, 2224. 

  Legendre, P. & Gallagher, E.D. (2001) Ecologically meaningful
  transformations for ordination of species data. \emph{Oecologia}
  \strong{129}, 271--280.

  Martino, C., Morton, J.T., Marotz, C.A., Thompson, L.R., Tripathi, A.,
  Knight, R. & Zengler, K. (2019) A novel sparse compositional technique
  reveals microbial perturbations.
  \emph{mSystems} \strong{4}, 1.

  Oksanen, J. (1983) Ordination of boreal heath-like vegetation with
  principal component analysis, correspondence analysis and
  multidimensional scaling. \emph{Vegetatio} \strong{52}, 181--189.
  }

\examples{
data(varespec)
sptrans <- decostand(varespec, "max")
apply(sptrans, 2, max)
sptrans <- wisconsin(varespec)

# CLR transformation for rows, with pseudocount
varespec.clr <- decostand(varespec, "clr", pseudocount=1)

# ALR transformation for rows, with pseudocount and reference sample
varespec.alr <- decostand(varespec, "alr", pseudocount=1, reference=1)

## Chi-square: PCA similar but not identical to CA.
## Use wcmdscale for weighted analysis and identical results.
sptrans <- decostand(varespec, "chi.square")
plot(procrustes(rda(sptrans), cca(varespec)))
}
\keyword{ multivariate}
\keyword{ manip }
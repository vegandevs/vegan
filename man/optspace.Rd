\encoding{UTF-8}
\name{optspace}
\alias{optspace}

\title{optspace: algorithm for matrix reconstruction from a partially revealed set}

\description{
This function was adapted from the original source code in the
\pkg{Roptspace} R package (version 0.2.3; MIT License) by Raghunandan
H. Keshavan, Andrea Montanari, Sewoong Oh (2010). See
\code{ROptSpace::OptSpace} for more information.  Let's assume an ideal
matrix \eqn{M} with \eqn{(m\times n)} entries with rank \eqn{r} and we
are given a partially observed matrix \eqn{M\_E} which contains many
missing entries. Matrix reconstruction - or completion - is the task
of filling in such entries. optspace is an efficient algorithm that
reconstructs \eqn{M} from \eqn{|E|=O(rn)} observed elements with
relative root mean square error (RMSE) \deqn{RMSE \le
C(\alpha)\sqrt{nr/|E|}}.
}

\usage{
optspace(x, ropt = 3, niter = 5, tol = 1e-5, verbose = FALSE)
}

\arguments{
  \item{x}{An \eqn{(n\times m)} matrix whose missing entries should be flagged as NA.}
  \item{ropt}{\code{FALSE} to guess the rank, or a positive integer as a pre-defined rank (default: 3).}
  \item{niter}{Maximum number of iterations allowed.}
  \item{tol}{Stopping criterion for reconstruction in Frobenius norm.}
  \item{verbose}{a logical value; \code{TRUE} to show progress, \code{FALSE} otherwise.}
}



\details{
This implementation removes the trimming step of the original
\code{Roptspace::OptSpace} code in order to leave feature filtering to the
user. Some of the defaults have been adjusted to better reflect
ecological data. The implementation has been adjusted for ecological
applications as in Martino et al. (2019). The imputed matrix (M) in
the optspace output includes matrix reconstruction (XSY'), with
subsequent centering for the columns and rows.
}

\value{Returns a named list containing:
  \item{X}{an \eqn{(n \times r)} matrix as left singular vectors.}
  \item{S}{an \eqn{(r \times r)} matrix as singular values.}
  \item{Y}{an \eqn{(m \times r)} matrix as right singular vectors.}    
  \item{dist}{a vector containing reconstruction errors at each successive iteration.}
  \item{M}{an \eqn{(n \times m)} imputed matrix, with columns and rows centered to zero.}        
}

\author{Leo Lahti and Cameron Martino, with adaptations of the method
implemented in \code{Roptspace::OptSpace} by Keshavan et
al. (2010).}

\references{ 

  Keshavan, R. H., Montanari, A., Oh, S. (2010).
  Matrix Completion From a Few Entries.
  \emph{IEEE Transactions on Information Theory} \strong{56}(6):2980--2998.

  Martino, C., Morton, J.T., Marotz, C.A., Thompson, L.R., Tripathi, A.,
  Knight, R. & Zengler, K. (2019) A novel sparse compositional technique
  reveals microbial perturbations.
  \emph{mSystems} \strong{4}, 1.

  }

\examples{

data(varespec)
# rclr transformation with no matrix completion for the 0/NA entries
x <- decostand(varespec, method = "rclr", impute = FALSE)
# Add matrix completion
xc <- optspace(x, ropt = 3, niter = 5, tol = 1e-5, verbose = FALSE)$M
}
\keyword{ multivariate}
\keyword{ manip }
\name{nullmodel}
\alias{nullmodel}
\alias{print.nullmodel}
\alias{simulate.nullmodel}
\alias{update.nullmodel}
\alias{str.nullmodel}
\alias{print.simmat}
\title{
Null Model and Simulation
}
\description{
The \code{nullmodel} function creates an object,
which can serve as a basis for Null Model simulation
via the \code{\link{simulate}} method.
The \code{\link{update}} method updates the nullmodel
object without sampling (effective for sequential algorithms).
}
\usage{
nullmodel(x, method)
\method{print}{nullmodel}(x, ...)
\method{simulate}{nullmodel}(object, nsim = 1, 
seed = NULL, burnin = 0, thin = 1, ...)
\method{update}{nullmodel}(object, nsim = 1, 
seed = NULL, ...)
\method{print}{simmat}(x, ...)
}
\arguments{
  \item{x}{
A community matrix.
For the \code{print} method, it is an object to be printed.
}
  \item{method}{
Character, specifying one of the null model algorithms
listed on the help page of \code{\link{commsim}}.
It can be a user supplied object of class \code{commsim}.
}
  \item{object}{
An object of class \code{nullmodel} returned by
the function \code{nullmodel}.
}
  \item{nsim}{
Positive integer, the number of simulated matrices to return.
For the \code{update} method, it is the number of
burnin steps made for sequential algorithms
to update the status of the input model \code{object}.
}
  \item{seed}{
An object specifying if and how the random number 
generator should be initialized ("seeded").
Either \code{NULL} or an integer that will be 
used in a call to \code{\link{set.seed}} before 
simulating the matrices. 
If set, the value is saved as the 
\code{"seed"} attribute of the returned value. 
The default, \code{NULL} will not change the 
random generator state, and return 
\code{\link{.Random.seed}} as the \code{"seed"}
 attribute, see Value. 
}
  \item{burnin}{
Nonnegative integer, specifying the number of steps
discarded before starting simulation.
Active only for sequential null model algorithms.
Ignored for non-sequential null model algorithms.
}
  \item{thin}{
Positive integer, number of simulation steps
made between each returned matrix.
Active only for sequential null model algorithms.
Ignored for non-sequential null model algorithms.
}
  \item{\dots}{
Additional arguments supplied to algorithms.
}
}
\details{
The purpose of the \code{nullmodel} function is to
create an object, where all necessary statistics of the
input matrix are calculated only once.
This information is reused, but not recalculated
in each step of the simulation process done by
the \code{simulate} method.

The \code{simulate} method carries out the simulation,
the simulated matrices are stored in an array.
For sequential algorithms, the method updates the state
of the input \code{nullmodel} object.
Therefore, it is possible to do diagnostic
tests on the returned \code{simmat} object,
and make further simulations, or use
increased thinning value if desired.

The \code{update} method makes burnin steps in case
of sequential algorithms to update the status of the 
input model without any attempt to return matrices.
For non-sequential algorithms the method does nothing.

\code{update} is the preferred way of making burnin iterations
without sampling. Alternatively, burnin can be done
via the \code{simulate} method. For convergence
diagnostics, it is recommended to use the
\code{simulate} method without burnin.
The input nullmodel object is updated, so further
samples can be simulated if desired without having
to start the process all over again. See Examples.
}
\value{
The function \code{nullmodel} returns an object of class \code{nullmodel}.
It is a set of objects sharing the same environment:
\itemize{
  \item{\code{data}: }{original matrix in integer mode.}
  \item{\code{nrow}: }{number of rows.}
  \item{\code{ncol}: }{number of columns.}
  \item{\code{rowSums}: }{row sums.}
  \item{\code{colSums}: }{column sums.}
  \item{\code{rowFreq}: }{row frequencies (number of nonzero cells).}
  \item{\code{colFreq}: }{column frequencies (number of nonzero cells).}
  \item{\code{totalSum}: }{total sum.}
  \item{\code{fill}: }{number of nonzero cells in the matrix.}
  \item{\code{commsim}: }{the \code{commsim} object as a result
    of the \code{method} argument.}
  \item{\code{state}: }{current state of the permutations, 
    a matrix similar to the original.
    It is \code{NULL} for non-sequential algorithms.}
  \item{\code{iter}: }{current number of iterations 
  for sequential algorithms.
    It is \code{NULL} for non-sequential algorithms.}
}

The \code{simulate} method returns an object of class \code{simmat}.
It is an array of simulated matrices (third dimension
corresponding to \code{nsim} argument).

The \code{update} method returns the current state (last updated matrix)
invisibly, and update the input object for sequential algorithms.
For non sequential algorithms, it returns \code{NULL}.
}
\author{
Jari Oksanen and Peter Solymos
}
\seealso{
\code{\link{commsim}}, \code{\link{make.commsim}}, 
\code{\link{permatfull}}, \code{\link{permatswap}}
}
\examples{
x <- matrix(rbinom(12*10, 1, 0.5)*rpois(12*10, 3), 12, 10)

## non-sequential nullmodel
(nm <- nullmodel(x, "r00"))
(sm <- simulate(nm, nsim=10))

## sequential nullmodel
(nm <- nullmodel(x, "swap"))
(sm1 <- simulate(nm, nsim=10, thin=5))
(sm2 <- simulate(nm, nsim=10, thin=5))

## sequential nullmodel with burnin and extra updating
(nm <- nullmodel(x, "swap"))
(sm1 <- simulate(nm, burnin=10, nsim=10, thin=5))
(sm2 <- simulate(nm, nsim=10, thin=5))

## sequential nullmodel with separate initial burnin
(nm <- nullmodel(x, "swap"))
nm <- update(nm, nsim=10)
(sm2 <- simulate(nm, nsim=10, thin=5))
}
\keyword{ multivariate }
\keyword{ datagen }

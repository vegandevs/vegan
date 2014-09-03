\name{RV}
\alias{RV}

\title{
RV coefficient
}
\description{
Calculates an RV coefficient.
}
\usage{
RV(X, Y)
}

\arguments{
  \item{X}{A matrix with the same number of rows as \code{Y}}
  \item{Y}{A matrix with the same number of rows as \code{X}}
}
\details{
 \code{RV} calculates the RV coefficient between two matrices.The RV coefficient is a multivariate generalization of the squared Pearson's correlation. In package \pkg{FactoMineR}, the function \code{\link[FactoMineR]{coeffRV}} calculates other variations of the RV coefficient and propose a test for it.  This function was designed to calculate only the basic RV coefficient making it much faster than \code{\link[FactoMineR]{coeffRV}}.
}
\value{
The RV coefficient. A numeric value ranging between 0 and 1
}
\references{
Escouffier, Y. (1973) Le traitement des variables vectorielles. Biometrics \strong{29}:751--760.\cr
}
\author{
F. Guillaume Blanchet
}

\seealso{
\code{\link[FactoMineR]{coeffRV}}
}
\examples{
data(beetle)
data(beetle.expl)

RV(beetle,beetle.expl)

}
\keyword{ multivariate }

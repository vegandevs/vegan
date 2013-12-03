\name{ccanova}
\alias{ccanova}
\alias{anova.prc}

\title{Permutation Test for Constrained Correspondence Analysis,
  Redundancy Analysis and Constrained Analysis of Principal Coordinates }
\description{
  This function is deprecated.
}
\usage{
ccanova(object, alpha=0.05, beta=0.01, step=100, perm.max=9999,
      by = NULL, ...)
}

\arguments{
  \item{object}{A result object from \code{\link{cca}}. }
  \item{alpha}{Targeted Type I error rate. }
  \item{beta}{Accepted Type II error rate. }
  \item{step}{Number of permutations during one step. }
  \item{perm.max}{Maximum number of permutations. }
  \item{by}{Setting \code{by = "axis"} will assess significance for each
    constrained axis, and setting \code{by = "terms"} will assess
    significance for each term (sequentially from first to last), and
    setting \code{by = "margin"} will assess the marginal effects of the
    terms (each marginal term analysed in a model with all other
    variables).}
  \item{\dots}{Parameters passed to other functions. 
    \code{anova.cca} passes all arguments to
    \code{permutest.cca}. In \code{anova} with \code{by = "axis"} you
    can use argument \code{cutoff} (defaults \code{1}) which stops
    permutations after exceeding the given level. }
}
\details{
  
  Function \code{ccanova} is a deprecated version of \code{anova.cca}.
  The new \code{link{anova.cca}} function is completely
  rewritten. With the same random number seed and the same number of
  permutations, the results are mostly identical. The only difference
  is that marginal tests (\code{by = "margin"}) are implemented
  differently. The most important difference is that the new code is
  based on the \pkg{permute} package and uses the fixed number of
  permutations in all tests (cf. below). Both the new and old
  functions are based \code{\link{permutest.cca}}, and the function
  can pass extra arguments to \code{permutest.cca}.

  In \code{ccanova} the number of permutations is controlled by
  targeted \dQuote{critical} \eqn{P} value (\code{alpha}) and accepted
  Type II or rejection error (\code{beta}).  If the results of
  permutations differ from the targeted \code{alpha} at risk level
  given by \code{beta}, the permutations are terminated.  If the
  current estimate of \eqn{P} does not differ significantly from
  \code{alpha} of the alternative hypothesis, the permutations are
  continued with \code{step} new permutations (at the first step, the
  number of permutations is \code{step - 1}).  However, with
  \code{by="terms"} a fixed number of permutations will be used, and
  this is given by argument \code{permutations}, or if this is
  missing, by \code{step}.
  
}

\note{
  The function is deprecated and replaced with \code{\link{anova.cca}}. 
}
\references{
  Legendre, P. and Legendre, L. (2012). \emph{Numerical Ecology}. 3rd
  English ed. Elsevier.
}
\author{Jari Oksanen}

\keyword{ multivariate }
\keyword{ htest }


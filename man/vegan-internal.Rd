\name{vegan-internal}
\alias{ordiGetData}
\alias{ordiParseFormula}
\alias{ordiNAexclude}
\alias{ordiNApredict}
\alias{getPermuteMatrix}
\alias{howHead}
\alias{centroids.cca}
\alias{ordiTerminfo}
\alias{pasteCall}
\alias{ordiArgAbsorber}
\alias{veganCovEllipse}
\alias{hierParseFormula}
\alias{veganMahatrans}

\title{Internal vegan functions}

\description{
  Internal vegan functions that are not intended to be called directly,
  but only within other functions.
}
\usage{
ordiGetData(call, env)
ordiParseFormula(formula, data, xlev = NULL,  na.action = na.fail,
    subset = NULL)
ordiTerminfo(d, data)
ordiNAexclude(x, excluded)
ordiNApredict(omit, x)
ordiArgAbsorber(..., shrink, origin, scaling, triangular,
                display, choices, const, FUN)
centroids.cca(x, mf, wt)
getPermuteMatrix(perm, N, strata = NULL)
howHead(x, ...)
pasteCall(call, prefix = "Call:")
veganCovEllipse(cov, center = c(0, 0), scale = 1, npoints = 100)
veganMahatrans(x, s2, tol = 1e-8)
hierParseFormula(formula, data)
}

\details{ The description here is only intended for \pkg{vegan}
  developers: these functions are not intended for users, but they
  only should be used within functions.  In general, these functions
  are not exported to the namespace, but you must use
  \code{\link{get}} or \code{\link{:::}} to directly call these
  functions.
  
  \code{ordiGetData} finds the model frame of constraints and
  conditions in constrained ordination in the defined
  \code{env}ironment.  \code{ordiParseFormula} returns a list of three
  matrices (dependent variables, and \code{\link{model.matrix}} of
  constraints and conditions, possibly \code{NULL}) needed in
  constrained ordination. Argument \code{xlev} is passed to
  \code{\link{model.frame}}. \code{ordiTermInfo} finds the term
  information for constrained ordination as described in
  \code{\link{cca.object}}. \code{ordiNAexclude} implements
  \code{na.action = na.exclude} for constrained ordination finding WA
  scores of CCA components and site scores of unconstrained component
  from \code{excluded} rows of observations. Function
  \code{ordiNApredict} pads the result object with these or with
  WA scores similarly as \code{\link{napredict}}.

  \code{ordiArgAbsorber} absorbs arguments of \code{\link{scores}}
  function of \pkg{vegan} so that these do not cause superfluous
  warnings in graphical function \code{FUN}. If you implement
  \code{scores} functions with new arguments, you should update
  \code{ordiArgAbsorber}. 

  \code{centroids.cca} finds the weighted centroids of variables.

  \code{getPermuteMatrix} interprets user input and returns a
  permutation matrix where each row gives indices of observations for
  a permutation. The input \code{perm} can be a single number for the
  number of simple permutations, a result of
  \code{\link[permute]{how}} defining a permutation scheme or a
  permutation matrix. Function \code{permuted.index} was used earlier
  to generate permutations, but it is now deprecated.

  \code{howHead} formats the permutation scheme of
  \code{\link[permute]{how}} for display. The formatting is more
  compact than the one used in \code{print} in the \pkg{permute}
  package, and shows only non-default choices. This output is normally
  used when printing the results of \pkg{vegan} permutations.

  \code{pasteCall} prints the function call so that it is nicely wrapped
  in \code{\link[utils]{Sweave}} output.

  \code{veganCovEllipse} finds the coordinates for drawing a
  covariance ellipse.

  \code{veganMahatrans} transforms data matrix so that its Euclidean
  distances are Mahalanobis distances. The input data \code{x} must be
  a matrix centred by columns, and \code{s2} its covariance matrix. If
  \code{s2} is not given, covariance matrix is found from \code{x}
  within the function.

  \code{hierParseFormula} returns a list of one matrix (left hand side)
  and a model frame with factors representing hierarchy levels 
  (right hand side) to be used in \code{\link{adipart}}, 
  \code{\link{multipart}} and \code{\link{hiersimu}}.
}

\keyword{internal }

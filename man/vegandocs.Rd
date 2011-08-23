\name{vegandocs}
\alias{vegandocs}

\title{Display Package Documentation}

\description{ Display NEWS, vignettes, other special documents or
  ChangeLog in \pkg{vegan}, or vignettes in \pkg{permute}. }

\usage{
vegandocs(doc = c("NEWS", "ChangeLog", "FAQ-vegan.pdf",
    "intro-vegan.pdf", "diversity-vegan.pdf", "decision-vegan.pdf",
    "partitioning.pdf", "permutations.pdf"))
}

\arguments{
  \item{doc}{The name of the document (partial match, case sensitive). }
}
\author{Jari Oksanen }

\note{ The function is a kluge, since \R does not have this facility
  (I hope it will come there). Function \code{\link{vignette}} only
  works with vignettes processed by \R, but the current function also
  shows other \code{pdf} documents. The \code{permutations.pdf}
  document is in the \pkg{permute} package.  }

\seealso{\code{\link{vignette}}. }
\examples{
\dontrun{
vegandocs("Change")
}
}
\keyword{ utilities }
\keyword{ documentation }

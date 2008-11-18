\name{vegandocs}
\alias{vegandocs}

\title{Display Package Documentation}
\description{
  Display package documentation using \code{pager} or \code{pdfviewer}
  defined in \code{\link{options}}.
}
\usage{
vegandocs(doc = c("NEWS", "ChangeLog", "FAQ-vegan.pdf",
    "intro-vegan.pdf", "diversity-vegan.pdf", "decision-vegan.pdf",
    "partitioning.pdf"))
}

\arguments{
  \item{doc}{The name of the document (partial match, case sensitive). }
}
\author{Jari Oksanen }
\note{
  The function is a kluge, since \R does not have this facility (I hope
  it will come there). Function \code{\link{vignette}} only works with
  vignettes. 
}
\seealso{\code{\link{vignette}}. }
\examples{
\dontrun{
vegandocs("Change")
}
}
\keyword{ utilities }
\keyword{ documentation }

\name{vegandocs}
\alias{vegandocs}

\title{Display Package Documentation}

\description{ Display NEWS, vignettes, other special documents or
  ChangeLog in \pkg{vegan}, or vignettes in \pkg{permute}. }

\usage{
vegandocs(doc = c("NEWS", "ONEWS", "ChangeLog", "FAQ-vegan.pdf",
    "intro-vegan.pdf", "diversity-vegan.pdf", "decision-vegan.pdf",
    "partitioning.pdf", "permutations.pdf"))
}

\arguments{
  \item{doc}{The name of the document (partial match, case sensitive). }
}

\details{
  You can read the following documents with this function:

  \itemize{
 
   \item \code{NEWS}: most important new functions, features, fixes
     etc. from the user's point of view. These can be also read using
     \R command \code{\link{news}} as \code{news(package = "vegan")}.

   \item \code{ONEWS}: old news about \pkg{vegan} version \code{1.*}
     before September 2011.

   \item \code{ChangeLog}: similar to news, but intended for
     developers wit more fine grained comments on internal changes.

   \item \code{FAQ-vegan}: Frequently Asked Questions. Consult here
     before writing to Mail groups.

   \item \code{intro-vegan}: a \code{\link{vignette}} demonstrating a
    simple, standard ordination analysis. This can be also read using
    \code{vignette("intro-vegan", package="vegan")}. 

   \item \code{diversity-vegan}: a \code{\link{vignette}} describing
     (most) diversity analyses in \pkg{vegan}. This can be also read
     using \code{vignette("diversity-vegan", package="vegan")}.

   \item \code{decision-vegan}: a \code{\link{vignette}} discussing
     design decisions in \pkg{vegan}.  Currently this discusses
     parallel processing in \pkg{vegan}, implementing nestedness
     temperature (\code{\link{nestedtemp}}), backtracking algorithm in
     community null models (\code{\link{make.commsim}}), scaling of
     RDA results, and why WA scores are used as default instead of LC
     scores in constrained ordination.

   \item \code{partitioning}: Detailed description of variation
     partitioning schemes used in \code{\link{varpart}}.

   \item \code{permutations}: a \code{\link{vignette}} in the
     \pkg{permute} package giving an introduction to restricted
     permutation schemes. You can also read this using
     \code{vignette("permutations", package="permute")}.

  }

}


\author{Jari Oksanen }

\note{ Function \code{\link{vignette}} only works with vignettes
  processed by \R, but the current function also shows other
  \code{pdf} documents. You can extract \R code from
  \code{\link{vignette}}s, but not from other documents (see
  \code{Examples}).

  The \code{permutations.pdf} document is in the \pkg{permute}
  package.

}

\seealso{\code{\link{vignette}}, \code{\link{news}}. }
\examples{
\dontrun{
## Read NEWS
vegandocs()
## Alternatively (perhaps with different formatting)
news(package="vegan")
## Read a vignette
vegandocs("intro")
## with vignette()
vignette("intro-vegan", package="vegan")
## extract R code
vig <- vignette("intro-vegan", package="vegan")
edit(vig)
##
}
}
\keyword{ utilities }
\keyword{ documentation }

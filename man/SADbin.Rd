\name{SADbin}
\alias{SADbin}
\title{
Species-abundance distribution (SAD) bins
}
\description{
This function constructs bins for a species-abundance distribution. It was designed so that every variation of binnings can easily be constructed. 
}
\usage{
SADbin(data, method = c("log", "modlog", "modhalflog"), base = 2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
A vector of species abundance for the whole study area.
}
  \item{method}{
A character string defining which binning method to use.
}
  \item{base}{
The log base to use to perform the binning. 
}
}
\details{
All value in \code{data} are expected to be count data. An error message will be sent if the value in \code{data} are not integers.

The different method of binning are defined as such :

log (using a base = 2) : Bin 1 = 1, Bin 2 = 2, Bin 3 = 3-4, Bin 4 = 5-8, ...

modlog (using a base = 2) : Bin 1 = 1, Bin 2 = 2-3, Bin 3 = 4-7, Bin 4 = 8-15,...

modhalflog (using a base = 2) : Bin 1 = half the number of species with 1 individuals, Bin 2  = half the number of species with 1 individuals and half of the number of species with 2 individuals, Bin 3 = half of the number of species with 2 individuals, all species with 3 individuals and half of the number of species with 4 individuals; and so on.
}

\value{
\item{bin}{A vector defining the number of species in each bin}
\item{sp.bin}{A matrix with the species in rows and the bins in columns. This matrix describe which species is in which bin}
}
\references{
Gray, J.S., A. Bjorgesaeter, and K. I. Ugland. 2006. On plotting species abundance distributions. \emph{Journal of Animal Ecology} \strong{75}:752--756.
}
\author{
F. Guillaume Blanchet
}

\examples{

data(mite)
SADmite<-SADbin(colSums(mite))
barplot(SADmite$bin,names.arg=paste("Bin",1:length(SADmite$bin)),las=2)

}
\keyword{ cluster }

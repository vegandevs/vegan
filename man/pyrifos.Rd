\name{pyrifos}
\alias{pyrifos} 
\docType{data}
\title{Response of Aquatic Invertebrates to Insecticide Treatment}
\description{
  The data are log transformed abundances of aquatic invertebrate in
  twelve ditches studied in eleven times before and after an insecticide
  treatment.
 
}
\usage{data(pyrifos)}
\format{
  A data frame with 132 observations on the log-transformed (\code{log(10*x + 1)}) abundances
  of 178 species. There are only twelve sites (ditches, mesocosms), but
  these were studied repeatedly in eleven occasions. The treatment
  levels, treatment times, or ditch ID's are not in the data frame, but
  the data are very regular, and the example below shows how to obtain
  these external variables.
}

\details{  This data set was obtained from an experiment in outdoor
  experimental ditches. Twelve mesocosms were allocated at random to
  treatments; four served as controls, and the remaining eight were 
  treated once with the insecticide chlorpyrifos, with nominal dose
  levels of 0.1, 0.9, 6, and 44  \eqn{\mu}{mu }g/ L in two mesocosms
  each. The example data set invertebrates.
  Sampling was done 11 times, from week -4 pre-treatment through
  week 24 post-treatment, giving a total of 132 samples (12 mesocosms
  times 11 sampling dates), see van den Brink & ter Braak (1999) for
  details. The data set contains only the species data,
  but the example below shows how to obtain the treatment, time and
  ditch ID variables. 
}
\source{
  CANOCO 4 example data, with the permission of Cajo J. F. ter Braak.
}
\references{
  van den Brink, P.J. & ter Braak, C.J.F. (1999). Principal response
  curves: Analysis of time-dependent multivariate responses of
  biological community to stress. Environmental Toxicology and
  Chemistry, 18, 138--148.
}
\examples{
data(pyrifos)
ditch <- gl(12, 1, length=132)
week <- gl(11, 12, labels=c(-4, -1, 0.1, 1, 2, 4, 8, 12, 15, 19, 24))
dose <- factor(rep(c(0.1, 0, 0, 0.9, 0, 44, 6, 0.1, 44, 0.9, 0, 6), 11))
}
\keyword{datasets}

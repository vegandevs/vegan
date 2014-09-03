\name{coeffCompare}
\alias{coeffCompare}

\title{Compare dissimilarities used within the RDA framework}
\description{
  This function compares association coefficients used through the RDA framework with a minimum spanning tree. It was designed to compare how information explained by one dissimilarity coefficient diverge from the information explained by another. The comparison is made simultaneously on the sites scores, the species scores and the canonical coefficients. 
}
\usage{
coeffCompare(ordires, ordisigniaxis, pval = 0.05)
}

\arguments{
  \item{ordires}{
 A list of \code{\link[vegan]{rda}} or \code{\link[vegan]{capscale}} result object that includes a series of RDA or distance-based RDA (db-RDA) performed with different association coefficients on the same data.
}
  \item{ordisigniaxis}{
 A list of \code{\link[vegan]{anova.cca}} object where each axis of each RDA (or db-RDA) given in \code{ordires} was tested for significance. This argument can also be a vector defining the number of significant axes in each RDA. See details. 
}
  \item{pval}{
 Numeric. P-value threshold to select the number of axes to use. This argument is only active if a list of \code{\link[vegan]{anova.cca}} object is given for the argument \code{ordisigniaxis}, otherwise it is not considered. Default is 0.05.
}
}
\details{
For the argument \code{ordisigniaxis}, if a vector defining the number of significant axes is given, it is assumed that the significant axes are selected in sequential order from the first axis.
The comparison made here rely on the RV coefficient \code{\link{RV}}, a multivariate generalization of the squared Pearson's correlation where matrix with the same number of rows are compared.
\code{coeffCompare} should be used prior to using \code{\link{consensusRDA}} because it informs the user about the different association coefficients considered interesting to perform a consensus RDA. An association coefficient presenting results too different from the others should not be included in the consensus RDA, it should be considered apart or discarded.
}
\value{
\item{RVmat}{A resemblance matrix of RV coefficients calculated from the sites scores matrices of RDA for all pairs of association coefficients.}
\item{mst}{minimum spanning tree calculated on (1-siteRVmat).}
}
\author{
F. Guillaume Blanchet
}

\seealso{
\code{\link{RV}}, \code{\link{consensusRDA}}
}
\examples{
###################################################################
### This example reproduces Figure 7b of Blanchet et al. (in press)
###################################################################
data(beetle)
data(beetle.expl)

### Construct results object
ndis<-10
ordiRes<-vector("list",length=ndis)

#---------------------------------------------
### Perform the various constrained ordination
#---------------------------------------------
### RDA species profile
sp<-beetle/apply(beetle,1,sum)
ordiRes[[1]]<-rda(sp~.,data=beetle.expl)
	
### RDA chord
chord<-beetle/sqrt(apply(beetle^2,1,sum))
ordiRes[[2]]<-rda(chord~.,data=beetle.expl)
	
### RDA Hellinger
hell<-decostand(beetle,method="hellinger")
ordiRes[[3]]<-rda(hell~.,data=beetle.expl)
	
### RDA chi2
chisq<-decostand(beetle,method="chi.square")
ordiRes[[4]]<-rda(chisq~.,data=beetle.expl)

### db-RDA Bray-Curtis
bray<-sqrt(vegdist(beetle,method="bray"))
ordiRes[[5]]<-capscale(bray~.,data=beetle.expl,comm=beetle)
	
### db-RDA square-root Bray-Curtis
bray.sqrt<-sqrt(vegdist(beetle^0.5,method="bray"))
ordiRes[[6]]<-capscale(bray.sqrt~.,data=beetle.expl,comm=beetle^0.5)
	
### db-RDA fourth-root Bray-Curtis
bray.fort<-sqrt(vegdist(beetle^0.25,method="bray"))
ordiRes[[7]]<-capscale(bray.fort~.,data=beetle.expl,comm=beetle^0.25)
	
### db-RDA modified Gower log 2
beetleLog2<-decostand(beetle, "log",logbase=2)
mGowerLog2<-vegdist(beetleLog2, "altGower")
ordiRes[[8]]<-capscale(mGowerLog2~.,data=beetle.expl,comm=beetleLog2)
	
### db-RDA modified Gower log 5
beetleLog5<-decostand(beetle, "log",logbase=5)
mGowerLog5<-vegdist(beetleLog5, "altGower")
ordiRes[[9]]<-capscale(mGowerLog5~.,data=beetle.expl,comm=beetleLog5)
	
### db-RDA modified Gower log 10
beetleLog10<-decostand(beetle, "log",logbase=10)
mGowerLog10<-vegdist(beetleLog10, "altGower")
ordiRes[[10]]<-capscale(mGowerLog10~.,data=beetle.expl,comm=beetleLog10)

### Compare association coefficients
AssoComp<-coeffCompare(ordiRes,rep(7,ndis))

#---------------------------------------------
### Draw a graphic to visualize the comparison
#---------------------------------------------
### Name of association coefficient compared
name<-c("Species profiles","Chord","Hellinger","Chi2","Bray-Curtis",
"(Bray-Curtis)^0.5","(Bray-Curtis)^0.25","mGowerlog2","mGowerlog5","mGowerlog10")

plot(AssoComp$mst,type="t",labels=name,xlab="",ylab="",main="MST Sites scores")

}
\keyword{ multivariate }

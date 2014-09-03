coeffCompare <-
function(ordires, ordisigniaxis,pval=0.05){
###
### Compare a series of association coefficient calculated through an RDA
###
### Arguments :
###
### ordires : A list of "cca rda" object from the vegan package gathering a series of RDA performed with different association coefficients on the same data.
### ordisigniaxis : A list of anova.cca object where each axis was tested for each RDA in ordires or a vector defining the number of significant axes in the RDA. See details 
### pval : Numeric. P-value threshold to select the number of axes to use. This argument is only active if a list of anova.cca object is given for the argument ordisigniaxis, otherwise it is not considered.
###
###
### Details :
### 
### For the argument ordisigniaxis, if a vector of number of significant axes is given, for each RDA, it is assumed that the significant axes are selected in sequential order from the first axis.  
### 
### Value :
###
### RVmat : A matrix of RV coefficients calculated for all pairs of association coefficients
### mst : minimum spanning tree calculated on (1-RVmat)
###
### F. Guillaume Blanchet - February 2012. (Modified November 2012, June 2013)
################################################################################
	#----------------
	#CC# Basic object
	#----------------
	#CC# Number of sites
	nsites<-nrow(scores(ordires[[1]],display="sites"))
	#CC# Number of juges (association coefficients)
	njuges<-length(ordires)
	
	#------------------
	#CC# General checks
	#------------------
	#### Check if ordires contains only RDAs
	allRDA<-sapply(ordires,function(x) any(class(x)=="rda"))
	if(!all(allRDA)){
		stop("One or more canonical ordination in 'ordires' is not an RDA")
	}
	
	#### Check if ordires and ordisigniaxis have the same number of components
	if(length(ordisigniaxis)!=length(ordires)){
		stop("'ordires' is not the same length as 'ordisigniaxis'")
	}
	
	ordisigniClass<-sapply(ordisigniaxis,function(x) class(x))
	
	#### Check if the capscale objects have the right number of species
	anycapscale<-which(sapply(ordires,function(x) any(class(x)=="capscale")))
	if(length(anycapscale) > 0){
		nspcapscale<-numeric()
		counter<-1
		for(i in anycapscale){
			nspcapscale[counter]<-nrow(scores(ordires[[i]],display="sp"))
			counter<-counter+1
		}
		if(any(nspcapscale==nsites)){
			stop("One or more of the analysis performed with capscale() did not include a site by species community matrix")
		}
	}
	
	#### If ordisigniaxis is a list
	if(is.list(ordisigniaxis)){
		#### Check P-values
		if(pval < 0 | pval > 1){
			stop("'pval' must range between 0 and 1")
		}
	
		allsigni<-sapply(ordisigniaxis,function(x) any(class(x)=="anova.cca"))
		if(!all(allsigni)){
			stop("One or more canonical ordination test in 'ordisigniaxis' is not an 'anova.cca' object")
		}
		
		anovaname<-unique(unlist(strsplit(sapply(ordisigniaxis,function(x) rownames(x)[1]),"1")))
		if(!all(anovaname == "RDA" | anovaname == "CAP")){
			stop("anova.cca by axis should be either 'RDA' or 'CAP'")
		}
		
		#CC# Extract the number of signicant axes to use for each canonical ordinations
		ordisigniaxis<-sapply(ordisigniaxis,function(x) length(which(x[,5]<=pval)))
	#### If ordisigniaxis is a vector
	}else{
		if(is.vector(ordisigniaxis)){
			if(!is.numeric(ordisigniaxis)){
				stop("'ordisigniaxis' should be numeric")
			}
		}else{
			stop("'ordisigniaxis' should either be a list of anova.cca objects or a vector of number of significant axes")
		}
	}
	
	#### If there are no significant axes
	if(any(ordisigniaxis < 1)){
		stop("One or more analysis does not have any significant axis remove it and start again")
	}
	
	#-------------------------------------------------------------
	#CC# Extract all the significant axes of matrix Z in scaling 1
	#-------------------------------------------------------------
	Zsigni<-vector("list",length=njuges)
	
	for(i in 1:njuges){
		Zsigni[[i]]<-scores(ordires[[i]],display="lc",choices=1:ordisigniaxis[i],scaling=1)
	}
	
	#CC# Coefficient RVf between the different association ceofficients
	RVassoCoeff<-matrix(NA,nrow=njuges,ncol=njuges)
	
	for(i in 1:njuges){
		for(j in 1:njuges){
			RVassoCoeff[i,j]<-RV(Zsigni[[i]],Zsigni[[j]])
		}
	}
	
	#CC# Construct a minimum spanning tree
	mst<-spantree(as.dist(1-RVassoCoeff))
	
	#CC# results
	res<-list(RVmat=RVassoCoeff,mst=mst)
	class(res)<-"coeffCompare"
	return(res)
}

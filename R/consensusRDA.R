consensusRDA <-
function(ordires, ordisigniaxis,resp.var,expl.var,pval=0.05,scaling=2){
###
### Calculates a consensus results for a series of canonical ordinations performed with different association coefficients on the same data. This function is only available for ordinations performed with canonical redundancy analysis (RDA) or variant of RDA such as distance-based RDA and transformation-based RDA
###
### Arguments :
###
### ordires : A list of "cca rda" object from the vegan package gathering a series of RDA performed with different association coefficients on the same data.
### ordisigniaxis : A list of anova.cca object where each axis was tested for each RDA in ordires or a vector defining the number of significant axes in the RDA. See details. 
### resp.var : Matrix of response variables
### expl.var : Matrix of explanatory variables
### pval : Numeric. P-value threshold to select the number of axes to use. This argument is only active if a list of anova.cca object is given for the argument ordisigniaxis, otherwise it is not considered.
### scaling : Type of scaling used to project the results. Default is 1 (distance).
###
###
### Details
### 
### For the argument ordisigniaxis, if a vector of number of significant axes is given, for each RDA, it is assumed that the significant axes are selected in sequential order from the first axis.  
### 
### Although it is possible to apply a scaling 3 to the RDA (it is available in the vegan package), this scaling should only be used for canonical correspondence analysis (CCA), it does not make any sense to use in the RDA framework.
###
### F. Guillaume Blanchet - February 2012 (Modified November 2012, June 2013)
################################################################################
	#----------------
	#CC# Basic object
	#----------------
	#CC# Number of sites
	nsites<-nrow(resp.var)
	#CC# Number of species
	nsp<-ncol(resp.var)
	#CC# Number of juges (association coefficients)
	njuges<-length(ordires)
	
	#------------------
	#CC# General checks
	#------------------
	#### Check if there are completely collinear explanatory variables and count the number of non-collinear variables
	expl.var.tmp<-expl.var
	deter<-det(cor(expl.var))
	
	#CC# threshold
	tresh<-10^-8
	
	while(deter<= tresh){
		expl.var.tmp<-expl.var.tmp[,-1]	
		deter<-det(cor(expl.var.tmp))
	}
	nNoncolldesc<-ncol(expl.var.tmp)
	
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
	
	#### Make sure that X is a matrix 
	expl.var<-as.matrix(expl.var)
		
	#-------------------------------------------------------------
	#CC# Extract all the significant axes of matrix Z in scaling 1
	#-------------------------------------------------------------
	Zsigni<-vector("list",length=njuges)
	Zsignimat<-matrix(NA,ncol=0,nrow=nsites)
	
	for(i in 1:njuges){
		Zsigni[[i]]<-scores(ordires[[i]],display="lc",choices=1:ordisigniaxis[i],scaling=1)
		Zsignimat<-cbind(Zsignimat,Zsigni[[i]])
	}
	
	#-------------------------------------------------
	#CC# Perform an RDA between Zsignimat and expl.var
	#-------------------------------------------------
	Zrda<-rda(Zsignimat,expl.var)
	ZrdaEigen<-eigenvals(Zrda)[1:nNoncolldesc]
	naxes<-length(ZrdaEigen)
	
	#CC# Extract the Z consensus results from RDA
	Zconsensus<-scores(Zrda,choice=1:nNoncolldesc,display="lc",scaling=1)
	
	#CC# Extract the C consensus results from RDA
	Cconsensus<-scores(Zrda,choice=1:nNoncolldesc,display="bp",scaling=1)
	
	#### This is the procedure proposed by Legendre and Legendre (2012, Subsection 9.3.3)
	#### It is also the procedure proposed by Oksanen et al. in vegan and described in Blanchet et al. (submitted)
	Uconsensus<-t(scale(resp.var,scale=FALSE))%*%Zconsensus%*%diag(ZrdaEigen^(-0.5))/sqrt(nsites-1)
	
	#----------------------
	#CC# Ordination Scaling
	#----------------------
	if(scaling==2){
		#-------------
		#CC# Consensus
		#-------------
		#CC# Sites
		Zconsensus<-sweep(Zconsensus,2,sqrt(ZrdaEigen/sum(ZrdaEigen)),"/")
		
		#CC# Species
		Uconsensus<-sweep(Uconsensus,2,sqrt(ZrdaEigen/sum(ZrdaEigen)),"*")
		
		#### No scaling is performed on descriptors
	}
	
	res<-list(values=ZrdaEigen,siteConsensus=Zconsensus,spConsensus=Uconsensus,descConsensus=Cconsensus)
	
	return(res)
}

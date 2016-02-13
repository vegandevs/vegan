coeffCompare <-
function(ordires, ordisigniaxis=NULL,pval=0.05,mst=TRUE,...){
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
	
        #-------------------------------------------------------------
        #CC# Perform ANOVA by axis if ordisigniaxis is not a vector of
        #CC# number of axes to consider
        #-------------------------------------------------------------
        if(!(is.vector(ordisigniaxis) && length(ordisigniaxis)==length(ordires))){
        	if(is.null(ordisigniaxis)){
				ordisigniaxis<-lapply(ordires, anova, by = "axis", cutoff = pval,...)
        	}else{
				stop("'ordisigniaxis' needs to be a vector of number of axes to consider for each RDA or defined as 'NULL'")
        	}
        }

	#### Check if ordires and ordisigniaxis have the same number of
    #### components
	if(length(ordisigniaxis)!=length(ordires)){
		stop("'ordires' is not the same length as 'ordisigniaxis'")
	}
	
	ordisigniClass<-sapply(ordisigniaxis,function(x) class(x))
	
	#### Check if the capscale objects have the right number of species
	spCheck<-sapply(ordires, function(x) nrow(scores(x,display="sp")))
	if(!all(spCheck==spCheck[1])){
		stop("The number of species differ for one analysis or one or more of the analysis performed with capscale() did not include a site by species matrix")
	}
	
	#### If ordisigniaxis is a list
	if(is.list(ordisigniaxis)){
		#### Check P-values
		if(pval < 0 | pval > 1){
			stop("'pval' must range between 0 and 1")
		}
	
		allsigni<-sapply(ordisigniaxis,function(x) any(class(x)=="anova"))
		if(!all(allsigni)){
			stop("One or more canonical ordination test in 'ordisigniaxis' is not an 'anova.cca' object")
		}
		
		anovaname<-unique(unlist(strsplit(sapply(ordisigniaxis,function(x) rownames(x)[1]),"1")))
		if(!all(anovaname == "RDA" | anovaname == "CAP")){
			stop("anova.cca by axis should be either 'RDA' or 'CAP'")
		}
		
		#CC# Extract the number of signicant axes to use for each
                #CC# canonical ordinations
		ordisigniaxis<-sapply(ordisigniaxis,function(x) length(which(x[,4]<=pval)))
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

	#### Check if the right side of the equation is the same for all 
	#### ordination in ordires
	checkY<-sapply(ordires,model.matrix)
	if(!is.matrix(checkY)){
		stop("One or more analysis does not have the same number of explanatory variables")
	}
	
	if(!all(checkY[,1]==checkY)){
		stop("One or more analysis was carried out on different explanatory variables")
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
	if(mst){
		mst<-spantree(as.dist(1-RVassoCoeff))
		#CC# results
		res<-list(RVmat=RVassoCoeff,mst=mst)
	}else{
		#CC# results
		res<-RVassoCoeff
	}
	class(res)<-"coeffCompare"
	return(res)
}

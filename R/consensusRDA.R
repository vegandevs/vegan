consensusRDA <-
function(ordires, ordisigniaxis,X,Y,pval=0.05,scaling=2){
	#----------------
	#CC# Basic object
	#----------------
	#CC# Number of sites
	nsites<-nrow(X)
	#CC# Number of species
	nsp<-ncol(X)
	#CC# Number of juges (association coefficients)
	njuges<-length(ordires)
	
	#------------------
	#CC# General checks
	#------------------
	#### Check if there are completely collinear explanatory variables 
	#### and count the number of non-collinear variables
	Y.tmp<-Y
	deter<-det(cor(Y))
	
	#CC# threshold
	tresh<-10^-8
	
	while(deter<= tresh){
		Y.tmp<-Y.tmp[,-1]	
		deter<-det(cor(Y.tmp))
	}
	nNoncolldesc<-ncol(Y.tmp)
	
	#### Check if ordires contains only RDAs
	allRDA<-sapply(ordires,function(x) any(class(x)=="rda"))
	if(!all(allRDA)){
		stop("One or more canonical ordination in 'ordires' is not an RDA")
	}
	
	#### Check if ordires and ordisigniaxis have the same number of 
	#### components
	if(length(ordisigniaxis)!=length(ordires)){
		stop("'ordires' is not the same length as 'ordisigniaxis'")
	}
	
	ordisigniClass<-sapply(ordisigniaxis, class)
	
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

		allsigni<-sapply(ordisigniaxis,function(x) any(class(x)=="anova"))
		if(!all(allsigni)){
			stop("One or more canonical ordination test in 'ordisigniaxis' is not an 'anova' object")
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
	
	#### Make sure that X is a matrix 
	Y<-as.matrix(Y)
	
	#### Check that all RDAs in ordires were carried out on the same
	#### Y data
	checkCall<-sapply(ordires,function(x) class(x$call[[2]]))
	if(!all(checkCall=="formula") | !all(checkCall=="call")){
		stop("All ordinations in 'ordires' needs to be performed the same way, using 'formulas' or 'matrix' not a combination of both")
	}
	
	if(all(checkCall=="formula")){
		
	}
	
	#### Check if the right side of the equation is the same for all 
	#### ordination in ordires
	checkY<-sapply(ordiRes,model.matrix)
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
	Zsignimat<-matrix(NA,ncol=0,nrow=nsites)
	
	for(i in 1:njuges){
		Zsigni[[i]]<-scores(ordires[[i]],display="lc",choices=1:ordisigniaxis[i],scaling=1)
		Zsignimat<-cbind(Zsignimat,Zsigni[[i]])
	}
	
	#------------------------------------------
	#CC# Perform an RDA between Zsignimat and Y
	#------------------------------------------
	Zrda<-rda(Zsignimat,Y)
	ZrdaEigen<-eigenvals(Zrda)[1:nNoncolldesc]
	naxes<-length(ZrdaEigen)
	
	#CC# Extract the Z consensus results from RDA
	Zconsensus<-scores(Zrda,choice=1:nNoncolldesc,display="lc",scaling=1)
	
	#CC# Extract the C consensus results from RDA
	Cconsensus<-scores(Zrda,choice=1:nNoncolldesc,display="bp",scaling=1)
	
	#### This is the procedure proposed by Legendre and 
	#### Legendre (2012, Subsection 9.3.3). It is also the procedure 
	#### proposed by Oksanen et al. in vegan and described in 
	#### Blanchet et al. (2014)
	Uconsensus<-t(scale(X,scale=FALSE))%*%Zconsensus%*%diag(ZrdaEigen^(-0.5))/sqrt(nsites-1)
	
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

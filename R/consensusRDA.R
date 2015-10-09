consensusRDA <-
function(ordires, X,Y,ordisigniaxis=NULL,pval=0.05,scaling=2,...){
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
        Q <- qr(Y)
        Y <- Y[, Q$pivot[seq_len(Q$rank)], drop = FALSE]
	
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
	
	ordisigniClass<-sapply(ordisigniaxis, class)
	
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
	
	#### Make sure that Y is a matrix 
	Y<-as.matrix(Y)
	
	#### Check that all RDAs in ordires were carried out on the same
	#### Y data
	checkCall<-sapply(ordires,function(x) class(x$call[[2]]))
	if(!all(checkCall=="formula")){
		if(!all(checkCall=="call")){
			stop("All ordinations in 'ordires' needs to be performed the same way, using 'formulas' or 'matrix' not a combination of both")
                }
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
	Zsignimat<-matrix(NA,ncol=0,nrow=nsites)
	
	for(i in 1:njuges){
		Zsigni[[i]]<-scores(ordires[[i]],display="lc",choices=1:ordisigniaxis[i],scaling=1)
		Zsignimat<-cbind(Zsignimat,Zsigni[[i]])
	}
	
	#------------------------------------------
	#CC# Perform an RDA between Zsignimat and Y
	#------------------------------------------
	Zrda<-rda(Zsignimat,Y)
        ZrdaAxes<-grep("RDA",names(eigenvals(Zrda)))
        ZrdaEigen<-eigenvals(Zrda)[ZrdaAxes]
	naxes<-length(ZrdaEigen)
	
	#CC# Extract the Z consensus results from RDA
	Zconsensus<-scores(Zrda,choice=ZrdaAxes,display="lc",scaling=1)
	
	#CC# Extract the C consensus results from RDA
	Cconsensus<-scores(Zrda,choice=ZrdaAxes,display="bp",scaling=1)
	
	#### This is the procedure proposed by Legendre and
        #### Legendre (2012, Subsection 9.3.3)  It is also the procedure
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

        class(res)<-"consensusRDA"
	return(res)
}

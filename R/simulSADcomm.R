simulSADcomm <-
function(sp.abund,expl.var,expl.rand.sel=TRUE,nexpl.comb=2,binary=FALSE,fix.expl=NULL,nsite=50,weight=NULL,range.weight=c(0,2),sd.expl=FALSE,norm=c(0,1)){
### Description:
### 
### Function that simulates data tables which have the same species
### abundance distribution.
###
### Arguments:
###
### sp.abund : A vector defining the number of species in a bin.
###            See "Details" for more information.
### nsite : Numeric. Number of sites (rows) in the resulting matrix. See details.
### expl.var : Matrix. Explanatory variables related to the species.
### expl.rand.sel : Logical. Whether explanatory should be randomly selected to construct species or a fixed combination should be given. (Default is TRUE)
### nexpl.comb : Numeric. The number of explanatory variables that will be combined together to construct the environmental variables. Default is 2.
### binary : Logical. Whether the site-by-species matrix is an abundance (FALSE) or a presence/absence (TRUE). Default is FALSE
### fix.expl : Matrix. Defines which combination of explanatory variables should be used to construct species. This argument is only active when expl.rand.sel=FALSE. See Details for more information.
### weight : Vector. Regression coefficient used to give weight on each species. If NULL weights are random selected through a random samping of a uniform distribution with a range defined by range.weight. Default is NULL.
### range.weight : Vector of length 2 giving the minimum and the maximum of a uniform distribution. This will be used to weight each species use to construct an explanatory variable. Default is 0 and 2.
### sd.expl : Logical. Whether the standard deviation of the Normal error is a multiplier of the standard deviation of the deterministic portion of the newly created explanatory variable (TRUE) or the pure standard deviation (FALSE). Default is FALSE.
### norm : Vector of length 2 giving the mean and a multiplier of the standard deviation of the deterministic portion of a newly created explanatory variable. Default is mean 0 and 1 time the standard deviation of the new deterministic explanatory variable.
###
### Details :
### 
### The argument "sp.abund" defines the species-abundance distribution structure of the data following the binnings proposed by Gray et al. (2006). For example, if the vector is (40,20,30), it means that there will be 40 species with 1 individual, 20 with 2 or 3 individuals, and 30 with 4 to 7 individuals.
###
### The individuals are assigned to the sites according to the the set of exlanatory variables given in expl.var. It is possible that a site occur with 0 individuals. They will be included in the analysis and dealt with a posteriori. 
###
### The explanatory variables are randomly sampled (without replacement) when combining (adding) explanatory variables together. The number of explanatory variables must be a multiple of nexpl.comb.
###
### Error was included to a species by multiplying a weight to the explanatory variable used to construct the species and by adding a normally distributed error term to that same explanatory variable. An error term with a standard deviation equal to the standard deviation of the explanatory variable allows for the explanatory variable to explain roughly 50% of the species it constructed.
###
### fix.expl is a matrix that has as many rows as there are species and as many columns as nexpl.comb (number of explanatory variables to combine). The numbers in fix.expl are integers that refers to the columns of expl.var. When fix.expl is used, nexpl.comb becomes meaningless.
###
### If a presence-absence matrix is constructed (binary=TRUE), sp.abund should be constructed in such a way that no bin should include species with an abundance larger than the number of sites. If it is not the case, an error message is sent. Within, this constraint, if the maximum of the last bin (the one with the largest abundance) is larger than the number of site, it will be automatically changed to the number of sites-1.
###
### Value :
###
### site.sp : The site (rows) by species (column) generated.
### sel.expl : A vector presenting the order explanatory variables used to model which species. The order follows the order of the species.
###
### Reference : 
### Gray, J. S., A. Bjorgeaeter, and K. I. Ugland. 2006. On plotting species abundance distributions, Journal of Animal Ecology. 75:752-756.
### 
### 
### F. Guillaume Blanchet - September 2010, July 2011
################################################################################
	if(!is.vector(sp.abund)){
		stop("'sp.abund' is not a vector")
	}
	
	is.wholenumber <- function(x, tol = .Machine$double.eps^0.5)  abs(x - round(x)) < tol

	nexpl.var<-ncol(expl.var)
	
	if(expl.rand.sel){
		nexpl.var.new<-ncol(expl.var)/nexpl.comb
		if(!is.wholenumber(nexpl.var.new)){
			stop("'expl.var' is not a multiple of 'nexp.comb'")
		}
	}
	
	if(nrow(expl.var)!=nsite){
		stop("'expl.var' should have the same number of row as 'nsite'")
	}
	
	#CC# Find the minimum and the maximium number of individuals for each bin define in sp.abund
	nbins<-length(sp.abund)
	min.bin<-2^(0:(nbins-1))
	max.bin<-2^(1:nbins)-1
	
	if(binary){
		if(max(min.bin) > nsite){
			stop("'sp.abund' has species with abundance too large")
		}
		max.bin[which.max(max.bin)]<-nsite-1
		
	}
	
	#CC# Construct site by specie result matrix
	nsp<-sum(sp.abund)
	site.sp<-matrix(0,nsite,nsp)
	
	#CC# Construct matrix presenting the environmental variable selection
	if(expl.rand.sel){
		expl.var.new<-matrix(NA,ncol=nexpl.var.new,nrow=nsite)
		sel.expl<-sample(1:ncol(expl.var))
		
		first<-seq(1,nexpl.var,by=nexpl.comb)
		last<-seq(nexpl.comb,nexpl.var,by=nexpl.comb)
		
		for(i in 1:nexpl.var.new){
			expl.var.new[,i]<-rowSums(expl.var[,sel.expl[first[i]:last[i]]])
		}
		
		sel.expl.new<-sample(1:nexpl.var.new,nsp,replace=TRUE)
		sd.expl.var.new<-apply(expl.var.new,2,sd)
	}else{
		expl.var.new<-matrix(NA,ncol=nsp,nrow=nsite)
		
		for(i in 1:nsp){
			expl.var.new[,i]<-rowSums(expl.var[,fix.expl[i,]])
			sel.expl.new<-1:nsp
			if(sd.expl){
				sd.expl.var.new<-apply(expl.var.new,2,sd)
			}else{
				sd.expl.var.new<-rep(1,ncol(expl.var.new))
			}
		}
	}
	
	#CC# Fill up site by species matrix
	sp<-1
	for(i in 1:nbins){
		if(sp.abund[i]>0){
			for(j in 1:sp.abund[i]){
				for(k in sample(min.bin[i]:max.bin[i],1)){
					#CC# Add a weight and an error term to the selected environmental variable
					error<-rnorm(nsite,mean=norm[1],sd=sd.expl.var.new[sel.expl.new[sp]]*norm[2])
					if(is.null(weight)){
						weight.rnd<-runif(1,range.weight[1],range.weight[2])
						smpl.prob<-abs(expl.var.new[,sel.expl.new[sp]]*weight.rnd+error)
						#CC# Consider the sign of the regression coefficient
						if(weight.rnd>0){
							smpl.prob<-smpl.prob/sum(smpl.prob)
						}else{
							smpl.prob<-(1/smpl.prob)/sum(1/smpl.prob)
						}
					}else{
						smpl.prob<-abs(expl.var.new[,sel.expl.new[sp]]*weight[sp]+error)
						#CC# Consider the sign of the regression coefficient
						if(weight[sp]>0){
							smpl.prob<-smpl.prob/sum(smpl.prob)
						}else{
							smpl.prob<-(1/smpl.prob)/sum(1/smpl.prob)
						}
					}
					#CC# Build presence/absence data
					if(binary){
						smpl.site<-sample(nsite,k,replace=FALSE,prob=smpl.prob)
						site.sp[smpl.site,sp]<-site.sp[smpl.site,sp]+1
					}else{
						smpl.site<-sample(nsite,k,replace=TRUE,prob=smpl.prob)
						for(l in smpl.site){
							site.sp[l,sp]<-site.sp[l,sp]+1
						}
					}
				}
				sp<-sp+1
			}
		}
	}
	
	if(expl.rand.sel){
		res<-list(site.sp,sel.expl)
		names(res)<-c("site.sp","sel.expl")
	}else{
		res<-list(site.sp,fix.expl)
		names(res)<-c("site.sp","sel.expl")
	}
	return(res)
}

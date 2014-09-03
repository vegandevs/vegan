SADbin <-
function(data,method=c("log","modlog","modhalflog"),base=2){
### Function that does method 1 of binning from Gray et al. (2006)
###
### Arguments:
###
### data: species abundant data
### method: binning method following Gray et al. (2006). Either "log", "modlog", or "modhalflog"
### base: base of the log to be used.
###
### copyleft - Guillaume Blanchet - August 2008
##################################################################
	method<-match.arg(method)
	"%w/o%" <- function(x,y) x[!x %in% y] #--  x without y
	
	#CC# General check
	if(any(data!=round(data))){
		stop("'data' should only include integer values")
	}
	
	#CC# Find the number of class to be used
	spmax<-max(data)

	#CC# Number of bins
	nbin<-ceiling(log(spmax+1,base=base))+1
	
	#CC# Find the levels of abundance of species
	lev<-as.numeric(levels(as.factor(data)))
	nlev<-nlevels(as.factor(data))

	#CC# Number of species per level
	nsp.lev<-vector(length=nlev)
	data.lev<-vector(length=length(data))
	
	for(i in 1:nlev){
		search<-which(data==lev[i])
		nsp.lev[i]<-length(search)
		data.lev[search]<-lev[i]
	}
	
	#### Find the number of species in each bins (and which species goes in which bins)
	#CC# Starting the species which give an integer with a log_base
	sp.div<-base^(1:(nbin-1))
	sp.div<-c(1,sp.div)
	
	bin.mat<-matrix(0,ncol=nbin,nrow=length(data))
	
	bin.sp.div<-vector(length=nbin)
	for(i in 1:nbin){
		if(length(which(lev==sp.div[i]))!=0){
			bin.sp.div[i]<-nsp.lev[which(lev==sp.div[i])]
			bin.mat[which(data.lev==sp.div[i]),i]<-1
		}
	}

	if(method=="modhalflog"){
		bin.sp.div2<-bin.sp.div/2
		bin.sp.div.good<-c(0,bin.sp.div2)+c(bin.sp.div2,0)
		
		bin.mat2<-bin.mat/2
		bin.mat<-cbind(0,bin.mat2)+cbind(bin.mat2,0)
	}

	#CC# Than with all the other ones
	spnot.div<-data %w/o% sp.div
	spnot.div2<-which((data %in% sp.div)==FALSE)
	
	bin.spnot.div<-vector(mode="numeric",length=nbin)

	if(method=="log" | method=="modhalflog"){
		for(i in 1:length(spnot.div)){
			bin.sel<-which(sp.div>spnot.div[i])[1]
			bin.spnot.div[bin.sel]<-bin.spnot.div[bin.sel]+1
			
			bin.mat[spnot.div2[i],bin.sel]<-1
		}
	}
	else if(method=="modlog"){

		for(i in 1:length(spnot.div)){
			bin.sel<-which(sp.div>spnot.div[i])[1]-1
			bin.spnot.div[bin.sel]<-bin.spnot.div[bin.sel]+1
			
			bin.mat[spnot.div2[i],bin.sel]<-1
		}
	}else{
		stop("method should be 'log' or 'modlog'")
	}
	
	#CC# Construct the resulting binning
	if(method=="modhalflog"){
		bin<-bin.sp.div.good[-length(bin.sp.div.good)]+bin.spnot.div
		bin.mat<-bin.mat[,-length(bin.sp.div.good)]
	}
	else if(method=="log"){
		bin<-bin.sp.div+bin.spnot.div
	}
	else if(method=="modlog"){
		bin<-bin.sp.div+bin.spnot.div
		bin<-bin[-length(bin)]
		bin.mat<-bin.mat[,-ncol(bin.mat)]
		
	}else{
		stop("method should be 'log' or 'modlog'")
	}
	
	return(list(bin=bin,sp.bin=bin.mat))
}

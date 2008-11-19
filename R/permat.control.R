`permat.control` <-
function(ptype="full", mtype="count", method="quasiswap", fixedmar="both", shuffle="ind", reg=NULL, hab=NULL, burnin=10000, thin=1000)
{
list(ptype=ptype, mtype=mtype, method=method, fixedmar=fixedmar, shuffle=shuffle, reg=reg, hab=hab, burnin=burnin, thin=thin)
}


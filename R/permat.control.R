`permat.control` <-
function(ptype="full", mtype="count", method="quasiswap", fixedmar="both", shuffle="both", strata=NULL, burnin=0, thin=1)
{
list(ptype=ptype, mtype=mtype, method=method, fixedmar=fixedmar, shuffle=shuffle, strata=strata, burnin=burnin, thin=thin)
}


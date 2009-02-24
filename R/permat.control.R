`permat.control` <-
function(ptype="full", mtype="count", method="quasiswap", fixedmar="both", shuffle="both", strata=NULL, burnin=10000, thin=1000)
{
list(ptype=ptype, mtype=mtype, method=method, fixedmar=fixedmar, shuffle=shuffle, strata=strata, burnin=burnin, thin=thin)
}


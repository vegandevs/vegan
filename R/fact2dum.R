# This function looks for factors in a data frame, extracts them and converts them to
# dummy binary variables, and reconstructs the data frame with the numerical variables
# first (if any) and the dummy variables after them.
# The dummy variables receive the names of the corresponding factor levels.
#
#  Daniel Borcard, Universite de Montreal
#  January 2018

fact2dum <- function(X, verbose = FALSE)
{

# Extract factors
X.fac <- as.data.frame(Filter(is.factor, X))
nfac <- length(names(X.fac))

# If no factor return data frame as is
if(nfac == 0)
{
	if(verbose) cat("There is no factor in the data frame, which is returned as is.")
    X.comb <- X
}

else
 {
	if(verbose) cat("The original data frame contains", nfac, "factor(s).\n")

# Extract non-factors from original data frame
X.num <- Filter(Negate(is.factor), X)
if(verbose)
{
   if(ncol(X.num) == 0)
   {
	   cat("No quantitative variable has been found in the data frame, only factors\n")
   }
}

# Empty matrix of dummy variables
X.dummy <- matrix(nrow = nrow(X.fac), ncol = 0)

# Produce dummy variables for each factor
for(i in 1 : ncol(X.fac))
   {
    class(X.fac[,i]) <- "factor"  # To avoid polynom. contrasts on ordered factors
   	nlev <- length(levels(X.fac[,i]))
  	tmp <- as.matrix(model.matrix(~X.fac[, i])[, -1])
   	colnames(tmp) <- levels(X.fac[, i])[-1]
   	X.dummy <- cbind(X.dummy, tmp)
   }

if(length(X.num) > 0)
{
   X.comb <- data.frame(X.num, X.dummy)
}
else
   X.comb <- data.frame(X.dummy)

ncontr <- ncol(X.dummy)
if(verbose)  cat("The output data frame contains", ncontr, "dummy variables in place of the", nfac, "factors.")
 }

X.comb

}

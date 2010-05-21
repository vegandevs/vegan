`mantel.correlog` <-
    function(D.eco, D.geo=NULL, XY=NULL, n.class=0, break.pts=NULL,
             cutoff=TRUE, r.type="pearson", nperm=999, mult="holm",
             progressive=TRUE)
{ 
    r.type <- match.arg(r.type, c("pearson", "spearman", "kendall"))
    mult   <- match.arg(mult, c("sidak", p.adjust.methods))
    
    epsilon <- .Machine$double.eps
    D.eco <- as.matrix(D.eco)

    ## Geographic distance matrix
    if(!is.null(D.geo)) {
	if(!is.null(XY))
            stop("You provided both a geographic distance matrix and a list of site coordinates. Which one should the function use?")
	D.geo <- as.matrix(D.geo)
    } else {
	if(is.null(XY)) {
            stop("You did not provide a geographic distance matrix nor a list of site coordinates")
        } else {
            D.geo <- as.matrix(dist(XY))
        }
    }
    
    n <- nrow(D.geo)
    if(n != nrow(D.eco))
	stop("Numbers of objects in D.eco and D.geo are not equal")
    n.dist <- n*(n-1)/2
    vec.D <- as.vector(as.dist(D.geo))
    vec.DD <- as.vector(D.geo)
    
    ## Number of classes and breakpoints
    
    if(!is.null(break.pts)) { 
	## Use the list of break points
	if(n.class > 0)
            stop("You provided both a number of classes and a list of break points. Which one should the function use?")
	n.class = length(break.pts) - 1 
        
    } else {
        ## No breakpoints have been provided: equal-width classes
        if(n.class == 0) { 
            ## Use Sturges rule to determine the number of classes
            n.class <- ceiling(1 + log(n.dist, base=2))
        }
        ## Compute the breakpoints from n.class
        start.pt <- min(vec.D)
        end.pt <- max(vec.D)
        width <- (end.pt - start.pt)/n.class
        break.pts <- vector(length=n.class+1)
        break.pts[n.class+1] <- end.pt
        for(i in 1:n.class) 
            break.pts[i] <- start.pt + width*(i-1) 
    }
    
    half.cl <- n.class %/% 2
    
    ## Move the first breakpoint a little bit to the left
    break.pts[1] <- break.pts[1] - epsilon   
    
    ## Find the break points and the class indices
    class.ind <- break.pts[1:n.class] +
        (0.5*(break.pts[2:(n.class+1)]-break.pts[1:n.class]))
    
    ## Create the matrix of distance classes
    vec2 <- vector(length=n^2)
    for(i in 1:n^2) 
        vec2[i] <- min( which(break.pts >= vec.DD[i]) ) - 1 
    
    ## Start assembling the vectors of results
    class.index <- NA
    n.dist <- NA
    mantel.r <- NA
    mantel.p <- NA
    ## check.sums = matrix(NA,n.class,1)

    ## Create a model-matrix for each distance class, then compute a Mantel test
    for(k in 1:n.class) {
	class.index <- c(class.index, class.ind[k])
	vec3 <- rep(0, n*n)
	sel <- which(vec2 == k)
	vec3[sel] <- 1
	mat.D2 <- matrix(vec3,n,n)
	diag(mat.D2) <- 0
	n.dis <- sum(mat.D2)
	n.dist <- c(n.dist, n.dis)
	if(n.dis == 0) {
            mantel.r <- c(mantel.r, NA)
            mantel.p <- c(mantel.p, NA)
        } else {
            row.sums <- apply(mat.D2, 1, sum)
            ## check.sums[k,1] = length(which(row.sums == 0))
            if((cutoff==FALSE) ||
               !(cutoff==TRUE && k > half.cl && any(row.sums == 0))) {
                temp <- mantel(mat.D2, D.eco, method=r.type, permutations=nperm)
                mantel.r <- c(mantel.r, -temp$statistic)
                temp.p <- temp$signif
                if(temp$statistic >= 0) {
                    temp.p <- ((temp.p*nperm)+1)/(nperm+1)
                } else {
                    temp.p <- (((1-temp.p)*nperm)+1)/(nperm+1)
                }
                mantel.p <- c(mantel.p, temp.p)
            } else {
                mantel.r <- c(mantel.r, NA)
                mantel.p <- c(mantel.p, NA)
            }
        }
    }
    
    mantel.res <- cbind(class.index, n.dist, mantel.r, mantel.p)
    mantel.res <- mantel.res[-1,]
    
    ## Note: vector 'mantel.p' starts with a NA value
    mantel.p <- mantel.p[-1]
    n.tests <- length(which(!is.na(mantel.p)))
    
    if(mult == "none") {
	colnames(mantel.res) <-
            c("class.index", "n.dist", "Mantel.cor", "Pr(Mantel)")
    } else {	
	## Correct P-values for multiple testing
        if(progressive) {
            p.corr <- mantel.p[1]
            if(mult == "sidak") {
                for(j in 2:n.tests) 
                    p.corr <- c(p.corr, 1-(1-mantel.p[j])^j)
            } else {
                for(j in 2:n.tests) {
                    temp <- p.adjust(mantel.p[1:j], method=mult)
                    p.corr <- c(p.corr, temp[j])
                }
            }
        } else {
            ## Correct all p-values for 'n.tests' simultaneous tests
            if(mult == "sidak") {
                p.corr <- 1 - (1 - mantel.p[1:n.tests])^n.tests
            } else {
                p.corr <- p.adjust(mantel.p[1:n.tests], method=mult)
            }
        }
	temp <- c(p.corr, rep(NA,(n.class-n.tests)))
	mantel.res <- cbind(mantel.res, temp)
	colnames(mantel.res) <-
            c("class.index", "n.dist", "Mantel.cor", "Pr(Mantel)", "Pr(corrected)")
    }
    rownames(mantel.res) <-
        rownames(mantel.res,do.NULL = FALSE, prefix = "D.cl.")
    
    ## Output the results
    res <- list(mantel.res=mantel.res, n.class=n.class, break.pts=break.pts,
                mult=mult, n.tests=n.tests, call=match.call() )
    class(res) <- "mantel.correlog"
    res
}

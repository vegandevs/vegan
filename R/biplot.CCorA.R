`biplot.CCorA` <-
    function(x, plot.type="ov", xlabs, plot.axes = 1:2, int=0.5, col.Y="red", col.X="blue", cex=c(0.7,0.9), ...)
{
    #### Internal function
	larger.frame <- function(mat, percent=0.10)
	# Produce an object plot 10% larger than strictly necessary
	{
	range.mat <- apply(mat,2,range)
	z <- apply(range.mat, 2, function(x) x[2]-x[1])
	range.mat[1,] <- range.mat[1,]-z*percent
	range.mat[2,] <- range.mat[2,]+z*percent
	range.mat
	}
	####
	
	TYPE <- c("objects","variables","ov","biplots")
	type <- pmatch(plot.type, TYPE)
	if(is.na(type)) stop("Invalid plot.type")

    epsilon <- sqrt(.Machine$double.eps)
	if(length(which(x$Eigenvalues > epsilon)) == 1)
		stop("Plot of axes (", paste(plot.axes, collapse=","),
			") not drawn because the solution has a single dimension.")
	if(max(plot.axes) > length(which(x$Eigenvalues > epsilon)))
		stop("Plot of axes (", paste(plot.axes, collapse=","),
			") not drawn because the solution has fewer dimensions.")

	if (missing(xlabs))
		xlabs <- rownames(x$Cy)
	else if (!is.null(xlabs) && is.na(xlabs))
		xlabs <- rep(NA, nrow(x$Cy))
	else if (is.null(xlabs))
		xlabs <- 1:nrow(x$Cy)
	#
	lf.Y <- larger.frame(x$Cy[,plot.axes])
	lf.X <- larger.frame(x$Cx[,plot.axes])
	#
	# Four plot types are available
	if(type == 1) {  # Object plots
		# cat('plot.type = objects')
		par(mfrow=c(1,2), pty = "s")
		plot(lf.Y, asp=1, xlab=colnames(x$Cy[,plot.axes[1]]), ylab=colnames(x$Cy[,plot.axes[2]]), type="n")
		points(x$Cy[,plot.axes], col=col.Y)   # Solid dot: pch=19 
		text(x$Cy[,plot.axes],labels=xlabs, pos=3, col=col.Y, cex=cex[1])
		title(main = c("CCorA object plot","First data table (Y)"), line=2)
		#
		plot(lf.X, asp=1, xlab=colnames(x$Cy[,plot.axes[1]]), ylab=colnames(x$Cy[,plot.axes[2]]), type="n")
		points(x$Cx[,plot.axes], col=col.X)   # Solid dot: pch=19 
		text(x$Cx[,plot.axes],labels=xlabs, pos=3, col=col.X, cex=cex[1])
		title(main = c("CCorA object plot","Second data table (X)"), line=2)
		###
		###
		} else if(type == 2) {  # Variable plots
		# cat('plot.type = variables')
		par(mfrow=c(1,2), pty = "s")
		plot(x$corr.Y.Cy[,plot.axes], asp=1, xlim=c(-1,1), ylim=c(-1,1), xlab=colnames(x$Cy[,plot.axes[1]]), 
			ylab=colnames(x$Cy[,plot.axes[2]]), type="n")
		text(x$corr.Y.Cy[,plot.axes],labels=rownames(x$corr.Y.Cy), pos=3, col=col.Y, cex=cex[2])
		arrows(0,0,x$corr.Y.Cy[,plot.axes[1]],x$corr.Y.Cy[,plot.axes[2]], length=0.05, col=col.Y)
		abline(h=0, v=0)
		lines(cos(seq(0, 2*pi, l=100)), sin(seq(0, 2*pi, l=100)))
		lines(int * cos(seq(0, 2*pi, l=100)), int * sin(seq(0, 2*pi, l=100)))
		title(main = c("CCorA variable plot","First data table (Y)"), line=2)
		#
		plot(x$corr.X.Cx[,plot.axes], asp=1, xlim=c(-1,1), ylim=c(-1,1), xlab=colnames(x$Cy[,plot.axes[1]]), 
			ylab=colnames(x$Cy[,plot.axes[2]]), type="n")
		text(x$corr.X.Cx[,plot.axes],labels=rownames(x$corr.X.Cx), pos=3, col=col.X, cex=cex[2])
		arrows(0,0,x$corr.X.Cx[,plot.axes[1]],x$corr.X.Cx[,plot.axes[2]], length=0.05, col=col.X)
		abline(h=0, v=0)
	    lines(cos(seq(0, 2*pi, l=100)), sin(seq(0, 2*pi, l=100)))
	    lines(int * cos(seq(0, 2*pi, l=100)), int * sin(seq(0, 2*pi, l=100)))
		title(main = c("CCorA variable plot","Second data table (X)"), line=2)
		###
		###
		} else if(type == 3) {  # Object and variable plots
		# cat('plot.type = ov')
		# par(mfrow=c(2,2), mar=c(4.5,3.5,2,1))
        layout(matrix(c(1,2,3,4), ncol = 2, nrow = 2, 
            byrow = TRUE), widths = 1, heights = c(0.5,0.5))
        par(pty = "s", mar = c(4.5,3.5,2,1))
		#
		plot(lf.Y, asp=1, xlab=colnames(x$Cy[,plot.axes[1]]), ylab=colnames(x$Cy[,plot.axes[2]]), type="n")
		points(x$Cy[,plot.axes], col=col.Y)   # Solid dot: pch=19 
		text(x$Cy[,plot.axes],labels=xlabs, pos=3, col=col.Y, cex=cex[1])
		title(main = c("First data table (Y)"), line=1)
		#
		plot(lf.X, asp=1, xlab=colnames(x$Cy[,plot.axes[1]]), ylab=colnames(x$Cy[,plot.axes[2]]), type="n")
		points(x$Cx[,plot.axes], col=col.X)   # Solid dot: pch=19 
		text(x$Cx[,plot.axes],labels=xlabs, pos=3, col=col.X, cex=cex[1])
		title(main = c("Second data table (X)"), line=1)
		#
		plot(x$corr.Y.Cy[,plot.axes], asp=1, xlim=c(-1,1), ylim=c(-1,1), xlab=colnames(x$Cy[,plot.axes[1]]), 
			ylab=colnames(x$Cy[,plot.axes[2]]), type="n")
		text(x$corr.Y.Cy[,plot.axes],labels=rownames(x$corr.Y.Cy), pos=3, col=col.Y, cex=cex[2])
		arrows(0,0,x$corr.Y.Cy[,plot.axes[1]],x$corr.Y.Cy[,plot.axes[2]], length=0.05, col=col.Y)
		abline(h=0, v=0)
		lines(cos(seq(0, 2*pi, l=100)), sin(seq(0, 2*pi, l=100)))
		lines(int * cos(seq(0, 2*pi, l=100)), int * sin(seq(0, 2*pi, l=100)))
		#
		plot(x$corr.X.Cx[,plot.axes], asp=1, xlim=c(-1,1), ylim=c(-1,1), xlab=colnames(x$Cy[,plot.axes[1]]), 
			ylab=colnames(x$Cy[,plot.axes[2]]), type="n")
		text(x$corr.X.Cx[,plot.axes],labels=rownames(x$corr.X.Cx), pos=3, col=col.X, cex=cex[2])
		arrows(0,0,x$corr.X.Cx[,plot.axes[1]],x$corr.X.Cx[,plot.axes[2]], length=0.05, col=col.X)
		abline(h=0, v=0)
	    lines(cos(seq(0, 2*pi, l=100)), sin(seq(0, 2*pi, l=100)))
	    lines(int * cos(seq(0, 2*pi, l=100)), int * sin(seq(0, 2*pi, l=100)))
		###
		###
		} else if(type == 4) {  # Biplots
		# cat('plot.type = biplot')
		par(mfrow=c(1,2), pty = "s")
		biplot(x$Cy[,plot.axes], x$corr.Y.Cy[,plot.axes], col=c("black",col.Y), xlim=lf.Y[,1], ylim=lf.Y[,2], 
			xlabs = xlabs, arrow.len=0.05, cex=cex, ...)
		title(main = c("CCorA biplot","First data table (Y)"), line=4)
		#
		biplot(x$Cx[,plot.axes], x$corr.X.Cx[,plot.axes], col=c("black",col.X), xlim=lf.X[,1], ylim=lf.X[,2], 
			xlabs = xlabs, arrow.len=0.05, cex=cex, ...)
		title(main = c("CCorA biplot","Second data table (X)"), line=4)
		}
	invisible()
}
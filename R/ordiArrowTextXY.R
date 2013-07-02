### Location of the text at the point of the arrow. 'vect' are the
### coordinates of the arrow heads, and 'labels' are the text used to
### label these heads, '...' passes arguments (such as 'cex') to
### strwidth() and strheight().
`ordiArrowTextXY` <-
    function (vect, labels, ...) 
{
    w <- strwidth(labels, ...)
    h <- strheight(labels, ...)
    ## slope of arrows
    b <- vect[,2]/vect[,1]
    ## offset based on string dimensions
    off <- cbind(sign(vect[,1]) * (w/2 + h/4), 0.75 * h * sign(vect[,2]))
    ## move the centre of the string to the continuation of the arrow
    for(i in 1:nrow(vect)) {
        move <- off[i,2] / b[i]
        ## arrow points to the top/bottom of the text box
        if (is.finite(move) && abs(move) <= abs(off[i, 1]))
            off[i, 1] <- move
        else {
            ## arrow points to a side of the text box
            move <- b[i] * off[i,1]
            off[i, 2] <- move
        }  
    }
    off + vect
}

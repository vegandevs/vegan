##############################################################
##  COMPUTES BEALS SMOOTHING FOR ALL SPECIES IN TABLE        #
##  This is a more complete function than the previous one   #
##  in the vegan package. The parameter values that give the #
##  equivalence are 'beals(x, NA, x, 0, include=TRUE)'       #
##                                                           #
##  'x'     matrix to be replaced by beals values            #
##  'reference' matrix to be used as source for joint occurrences#
##  'type'  sets the way to use abundance values             #
##       0 - presence/absence                                #
##       1 - abundances for conditioned probabilities        #
##       2 - abundances for weighted average                 #
##       3 - abundances for both                             #
##  'species' a column index used to compute Beals function  #
##        for a single species. The default (NA) indicates   #
##        all species.                                       #
##  'include' flag to include target species in the computation#
##############################################################
beals<-function(x, species=NA, reference=x, type=0, include=TRUE)
{
    refX <- reference
    # this checks whether it was chosen from available options
    mode <- as.numeric(match.arg(as.character(type), c("0","1","2","3")))
    spIndex <- species
    incSp <- include
#    if(is.null(refX))              # I think this is unnecessary because of the default above
#            refX<-x                # and it should look like missing(refX) and not is.null(refX)
   refX <- as.matrix(refX)
   x <- as.matrix(x)
   if(mode==0 || mode ==2) refX <- ifelse(refX > 0, 1, 0)
   if(mode==0 || mode ==1) x <- ifelse(x > 0, 1, 0)
   #Computes conditioned probabilities
   if(is.na(spIndex)){
       M <- crossprod(ifelse(refX > 0, 1, 0),refX)
      C <-diag(M)
     M <- sweep(M, 2, replace(C,C==0,1), "/")
     if(!incSp)
         for (i in 1:ncol(refX))
            M[i,i]<-0
    } else {
       C<-colSums(refX)
       M<-crossprod(refX,ifelse(refX > 0, 1, 0)[,spIndex])
       M<-M/replace(C,C==0,1)
       if(!incSp)
          M[spIndex]<-0
    }
   #Average of conditioned probabilities
   S <- rowSums(x)
   if(is.na(spIndex)) {
       b <-x
       for (i in 1:nrow(x)) {
       b[i, ] <- rowSums(sweep(M, 2, x[i, ], "*"))
       }
       SM<-rep(S,ncol(x))
       if(!incSp)
           SM<-SM-x
       b <- b/replace(SM,SM==0,1)
   } else {
       b<-rowSums(sweep(x,2,M,"*"))
       if(!incSp)
           S<-S-x[,spIndex]
       b <- b/replace(S,S==0,1)
   }
#   attr(b, "smoothtype") <- mode
#   class(b) <- c("beals", class(b))        # this can later be used to write the methods, i.e. beals test, etc.
   return(b)
}

\name{NEWS}
\title{vegan News}
\encoding{UTF-8}

\section{Changes in version 2.0-beta}{

  \subsection{GENERAL}{
    \itemize{
      
      \item Peter Minchin joins the \pkg{vegan} team.
      
      \item \pkg{vegan} implements standard \R \code{NAMESPACE}. In
      general, \code{S3} methods are not exported which means that you
      cannot directly use or see contents of functions like
      \code{cca.default}, \code{plot.cca} or \code{anova.cca}. To use
      these functions you should rely on \R delegation and simply use
      \code{cca} and for its result objects use \code{plot} and
      \code{anova} without suffix \code{.cca}. To see the contents of
      the functión you can use \code{:::}, such as
      \code{vegan:::cca.default}. This change may break packages,
      documents or scripts that rely on non-exported names.
      
      \item \pkg{vegan} depends on the \pkg{permute} package. This
      package provides powerful tools for restricted permutation
      schemes. All \pkg{vegan} permutation will gradually move to use
      \pkg{permute}, but currently only \code{betadisper} uses the new
      feature. 
      
    }
  } % end GENERAL
  
  \subsection{NEW FUNCTIONS}{
    \itemize{
      
      \item \code{monoMDS}: a new function for non-metric
      multidimensional scaling (NMDS). This function replaces
      \code{MASS::isoMDS} as the default method in \code{metaMDS}. Major
      advantages of \code{monoMDS} is that it has so called
      \sQuote{primary} or \sQuote{weak} tie treatment which means that
      it can split tied observed dissimilarities.  \sQuote{Weak} tie
      treatment improves ordination of heterogeneous data sets, because
      maximum dissimilarities of \eqn{1} can be split. In addition to
      global NMDS, \code{monoMDS} can perform local and hybrid NMDS and
      metric MDS.  It can also handle missing and zero dissimilarities.
      Moreover, \code{monoMDS} is faster than previous alternatives.
      
      \item \code{fitspecaccum} to fit non-linear regression models to
      the species accumulation results from \code{specaccum}. The
      function can use new self-starting species accumulation models
      in \pkg{vegan} or other self-starting non-linear regression
      models in \R. The functin can fit Arrhenius, Gleason, Gitay,
      Lomolino (in \pkg{vegan}), asymptotic, Gompertz,
      Michaelis-Menten, logistic and Weibull (in base \R) models.
      
      \item Self-starting non-linear species accumulation models
      \code{SSarrhenius}, \code{SSgleason}, \code{SSgitay} and
      \code{SSlomolino}. These can be used with \code{fitspecaccum} or
      directly in non-linear regression with \code{nls}. 
      
    }
  } % end NEW FUNCTIONS
  
  \subsection{NEW FEATURES}{
    
    \itemize{
      
      \item \code{betadisper} uses restricted permutation of the
      \pkg{permute} package. 
      
      \item \code{metaMDS} uses \code{monoMDS} as its default ordination
      engine. Function gains new argument \code{engine} that can be used
      to alternatively select \code{MASS::isoMDS}.  The default is not
      to use \code{stepacross} with \code{monoMDS} because its
      \sQuote{weak} tie treatment can cope with tied maximum
      dissimilarities of one. However, \code{stepacross} is always used
      with \code{isoMDS} because it cannot handle adequately these tied
      maximum dissimilarities. 
      
    }
  } % end NEW FEATURES
  
  \subsection{BUG FIXES}{
    \itemize{
      \item None yet.
    }
  } % end BUG FIXES
  
  \subsection{ANALYSES}{
    \itemize{
      
      \item \code{treeheight} uses much snappier code. The results
      should be unchanged.
      
    }
  } % end ANALYSES
}% end VERSTION 2.0   
   	  

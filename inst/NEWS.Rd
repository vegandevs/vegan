\name{NEWS}
\title{vegan News}
\encoding{UTF-8}

\section{Changes in version 2.0-0}{
  \subsection{GENERAL}{
    \itemize{
      
      \item Peter Minchin joins the \pkg{vegan} team.
      
      \item \pkg{vegan} implements standard \R \code{NAMESPACE}. In
      general, \code{S3} methods are not exported which means that you
      cannot directly use or see contents of functions like
      \code{cca.default}, \code{plot.cca} or \code{anova.ccabyterm}. To
      use these functions you should rely on \R delegation and simply
      use \code{cca} and for its result objects use \code{plot} and
      \code{anova} without suffix \code{.cca}. To see the contents of
      the function you can use \code{:::}, such as
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
      \code{MASS::isoMDS} as the default method in
      \code{metaMDS}. Major advantages of \code{monoMDS} are that it
      has \sQuote{weak} (\sQuote{primary}) tie treatment which means
      that it can split tied observed dissimilarities.  \sQuote{Weak}
      tie treatment improves ordination of heterogeneous data sets,
      because maximum dissimilarities of \eqn{1} can be split. In
      addition to global NMDS, \code{monoMDS} can perform local and
      hybrid NMDS and metric MDS.  It can also handle missing and zero
      dissimilarities.  Moreover, \code{monoMDS} is faster than
      previous alternatives. The function uses \code{Fortran} code
      written by Peter Minchin.

      \item \code{MDSrotate} a new function to replace
      \code{metaMDSrotate}. This function can rotate both \code{metaMDS}
      and \code{monoMDS} results so that the first axis is parallel to
      an environmental vector.

      \item \code{eventstar} finds the minimum of the evenness profile
      on the Tsallis entropy, and uses this to find the corresponding
      values of diversity, evenness and numbers equivalent following
      Mendes et al. (\emph{Ecography} 31, 450-456; 2008). The code was
      contributed by Eduardo Ribeira Cunha and Heloisa Beatriz Antoniazi
      Evangelista and adapted to \pkg{vegan} by Peter Solymos.
      
      \item \code{fitspecaccum} fits non-linear regression models to
      the species accumulation results from \code{specaccum}. The
      function can use new self-starting species accumulation models
      in \pkg{vegan} or other self-starting non-linear regression
      models in \R. The function can fit Arrhenius, Gleason, Gitay,
      Lomolino (in \pkg{vegan}), asymptotic, Gompertz,
      Michaelis-Menten, logistic and Weibull (in base \R) models. The
      function has \code{plot} and \code{predict} methods.
      
      \item Self-starting non-linear species accumulation models
      \code{SSarrhenius}, \code{SSgleason}, \code{SSgitay} and
      \code{SSlomolino}. These can be used with \code{fitspecaccum} or
      directly in non-linear regression with \code{nls}. These functions
      were implemented because they were found good for species-area
      models by Dengler (\emph{J. Biogeogr.} 36, 728-744; 2009).
      
    }
  } % end NEW FUNCTIONS
  
  \subsection{NEW FEATURES}{
    \itemize{

      \item \code{adonis}, \code{anosim}, \code{meandist} and
      \code{mrpp} warn on negative dissimilarities, and
      \code{betadisper} refuses to analyse them. All these functions
      expect dissimilarities, and giving something else (like
      correlations) probably is a user error.
      
      \item \code{betadisper} uses restricted permutation of the
      \pkg{permute} package. 
      
      \item \code{metaMDS} uses \code{monoMDS} as its default ordination
      engine. Function gains new argument \code{engine} that can be used
      to alternatively select \code{MASS::isoMDS}.  The default is not
      to use \code{stepacross} with \code{monoMDS} because its
      \sQuote{weak} tie treatment can cope with tied maximum
      dissimilarities of one. However, \code{stepacross} is the default
      with \code{isoMDS} because it cannot handle adequately these tied
      maximum dissimilarities.
      
      \item \code{specaccum} gained \code{predict} method which uses
      either linear or spline interpolation for data between observed
      points. Extrapolation is possible with spline interpolation, but
      may make little sense.
      
      \item \code{specpool} can handle missing values or empty factor
      levels in the grouping factor \code{pool}.  Now also checks that
      the length of the \code{pool} matches the number of
      observations.

    }
  } % end NEW FEATURES

  \subsection{DEPRECATED AND DEFUNCT}{
    \itemize{

      \item \code{metaMDSrotate} was replaced with \code{MDSrotate}
      that can also handle the results of \code{monoMDS}.

      \item \code{permuted.index2} and other \dQuote{new} permutation
      code was removed in favour of the \pkg{permute} package. This code
      was not intended for normal use, but packages depending on that
      code in \pkg{vegan} should instead depend on \pkg{permute}.

    }  
    
  } % end DEPRECATED
  
  \subsection{ANALYSES}{
    \itemize{
      
      \item \code{treeheight} uses much snappier code. The results
      should be unchanged.
      
    }
  } % end ANALYSES
}% end VERSION 2.0   
   	  

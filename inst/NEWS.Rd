\name{NEWS}
\title{vegan News}
\encoding{UTF-8}

\section{Changes in version 2.0-3}{

  \subsection{NEW FUNCTIONS}{
    \itemize{

      \item Added new nestedness functions \code{nestedbetasor} and
      \code{nestedbetajac} that implement multiple-site dissimilarity
      indices and their decomposition into turnover and nestedness
      components following Baselga (\emph{Global Ecology and
      Biogeography} 19, 134--143; 2010).

      \item Added function \code{rarecurve} to draw rarefaction curves
      for each row (sampling unit) of the input data, optionally with
      lines showing rarefied species richness with given sample size
      for each curve.

      \item Added function \code{simper} that implements
      \dQuote{similarity percentages} of Clarke (\emph{Australian
      Journal of Ecology} 18, 117--143; 1993).  The method compares
      two or more groups and decomposes the average between-group
      Bray-Curtis dissimilarity index to contributions by individual
      species.  The code was developed in 
      \href{https://github.com/jarioksa/vegan}{GitHub}
      by Eduard Szöcz (Uni Landau, Germany).

    }
  } % end new functions

  \subsection{BUG FIXES}{
    \itemize{

      \item \code{betadisper()} failed when the \code{groups} was a
      factor with empty levels.

      \item Some constrained ordination methods and their support
      functions are more robust in border cases (completely aliased
      effects, saturated models, user requests for non-existng scores
      etc). Concerns \code{capscale}, \code{ordistep}, \code{varpart},
      \code{plot} function for constrained ordination, and
      \code{anova(<cca.object>, by = "margin")}.

      \item The \code{scores} function for \code{monoMDS} did not
      honour \code{choices} argument and hence dimensions could not be
      chosen in \code{plot}.

      \item The default \code{scores} method failed if the number of
      requested axes was higher than the ordination object had.  This
      was reported as an error in \code{ordiplot} in
      \href{https://stat.ethz.ch/pipermail/r-sig-ecology/2012-February/002768.html}{R-sig-ecology} mailing list.

    }
  } % end bug fixes	
		  
  \subsection{NEW FEATURES}{
     \itemize{

       \item \code{metaMDS} argument \code{noshare = 0} is now
       regarded as a numeric threshold that always triggers extended
       dissimilarities (\code{stepacross}), instead of being treated
       as synonymous with \code{noshare = FALSE} which always
       suppresses extended dissimilarities.

       \item Nestedness discrepancy index \code{nesteddisc} gained a
       new argument that allows user to set the number of iterations
       in optimizing the index.

       \item \code{oecosimu} displays the mean of simulations and
       describes alternative hypothesis more clearly in the printed
       output. 

       \item Implemented adjusted \eqn{R^2}{R-squared} for partial
       RDA. For partial model \code{rda(Y ~ X1 + Condition(X2))} this
       is the same as the component \code{[a] = X1|X2} in variance
       partition in \code{varpart} and describes the marginal (unique)
       effect of constraining term to adjusted \eqn{R^2}{R-squared}.

       \item Added Cao dissimilarity (CYd) as a new dissimilarity
       method in \code{vegdist} following Cao et al., \emph{Water
       Envir Res} 69, 95--106 (1997). The index should be good for
       data with high beta diversity and variable sampling
       intensity. Thanks to consultation to Yong Cao (Univ Illinois,
       USA).
     
     }
  } % end new features
} % end version 2.0-3

\section{Changes in version 2.0-2}{
  \subsection{BUG FIXES}{
   \itemize{

     \item Function \code{capscale} failed if constrained component
     had zero rank. This happened most likely in partial models when
     the conditions aliased constraints. The problem was observed in
     \code{anova(..., by ="margin")} which uses partial models to
     analyses the marginal effects, and was reported in an email
     message to
     \href{https://stat.ethz.ch/pipermail/r-help/2011-October/293077.html}{R-News
     mailing list}.

     \item \code{stressplot} and \code{goodness} sometimes failed when
     \code{metaMDS} was based on \code{isoMDS} (\pkg{MASS} package)
     because \code{metaMDSdist} did not use the same defaults for
     step-across (extended) dissimilarities as \code{metaMDS(...,
     engine = "isoMDS")}. The change of defaults can also influence
     triggering of step-across in \code{capscale(..., metaMDSdist =
     TRUE)}.

     \item \code{adonis} contained a minor bug resulting from
     incomplete implementation of a speed-up that did not affect the
     results. In fixing this bug, a further bug was identified in
     transposing the hat matrices. This second bug was only active
     following fixing of the first bug. In fixing both bugs, a
     speed-up in the internal f.test() function is fully
     realised. Reported by Nicholas Lewin-Koh.
   
     }
 }  % end bug fixes 
 
 \subsection{NEW FEATURES}{
   \itemize{

     \item \code{ordiarrows} and \code{ordisegments} gained argument
     \code{order.by} that gives a variable to sort points within
     \code{groups}. Earlier the points were assumed to be in order.
      
     \item Function \code{ordispider} invisibly returns the
     coordinates to which the points were connected.  Typically these
     are class centroids of each point, but for constrained ordination
     with no \code{groups} they are the LC scores.

    }		 
  }  %end new features
		  
} %end version 2.0-2

\section{Changes in version 2.0-1}{
  \subsection{NEW FUNCTIONS}{
    \itemize{
  
      \item \code{clamtest}: new function to classify species as
      generalists and specialists in two distinct habitats (CLAM test of
      Chazdon et al., \emph{Ecology} 92, 1332--1343; 2011).  The test is
      based on multinomial distribution of individuals in two habitat
      types or sampling units, and it is applicable only to count data
      with no over-dispersion.

      \item \code{as.preston} gained \code{plot} and \code{lines}
      methods, and \code{as.fisher} gained \code{plot} method (which
      also can add items to existing plots). These are similar as
      \code{plot} and \code{lines} for \code{prestonfit} and
      \code{fisherfit}, but display only data without the fitted lines.

      \item \code{raupcrick}: new function to implement Raup-Crick
      dissimilarity as a probability of number of co-occurring species
      with occurrence probabilities proportional to species
      frequencies.  \pkg{Vegan} has Raup-Crick index as a choice in
      \code{vegdist}, but that uses equal sampling probabilities for
      species and analytic equations. The new \code{raupcrick}
      function uses simulation with \code{oecosimu}. The function
      follows Chase et al. (2011) \emph{Ecosphere} 2:art24
      [\href{http://www.esajournals.org/doi/abs/10.1890/ES10-00117.1}{doi:10.1890/ES10-00117.1}],
      and was developed with the consultation of Brian Inouye.

    }
  } % end NEW FUNCTIONS

  \subsection{BUG FIXES}{
    \itemize{

      \item Function \code{meandist} could scramble items and give
      wrong results, especially when the \code{grouping} was
      numerical. The problem was reported by Dr Miguel Alvarez
      (Univ. Bonn).

      \item \code{metaMDS} did not reset \code{tries} when a new model
      was started with a \code{previous.best} solution from a different
      model. 

      \item Function \code{permatswap} for community null models using
      quantitative swap never swapped items in a \eqn{2 \times 2}{2 by 2}
      submatrix if all cells were filled.

      \item The result from \code{permutest.cca} could not be
      \code{update}d because of a \file{NAMESPACE} issue.

      \item \R 2.14.0 changed so that it does not accept using
      \code{sd()} function for matrices (which was the behaviour at
      least since \R 1.0-0), and several \pkg{vegan} functions were
      changed to adapt to this change (\code{rda}, \code{capscale},
      \code{simulate} methods for \code{rda}, \code{cca} and
      \code{capscale}). The change in \R 2.14.0 does not influence the
      results but you probably wish to upgrade \pkg{vegan} to avoid
      annoying warnings.

    }	
  } % end BUG FIXES

  \subsection{ANALYSES}{
    \itemize{

      \item \code{nesteddisc} is slacker and hence faster when trying
      to optimize the statistic for tied column frequencies. Tracing
      showed that in most cases an improved ordering was found rather
      early in tries, and the results are equally good in most cases.

    }
  } % end ANALYSES 

} % end version 2.0-1

\section{Changes in version 2.0-0}{
  \subsection{GENERAL}{
    \itemize{
      
      \item Peter Minchin joins the \pkg{vegan} team.
      
      \item \pkg{vegan} implements standard \R \file{NAMESPACE}. In
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
   	  

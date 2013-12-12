\name{NEWS}
\title{vegan News}
\encoding{UTF-8}

\section{Changes in version 2.0-10}{

  \subsection{GENERAL}{
    \itemize{

      \item This version is adapted to the changes in \pkg{permute}
      package version 0.8-0 and no more triggers NOTEs in package
      checks.  This release may be the last of the 2.0 series, and the
      next \pkg{vegan} release is scheduled to be a major release with
      newly designed \code{oecosimu} and community pattern simulation,
      support for parallel processing, and full support of the
      \pkg{permute} package. If you are interested in these
      developments, you may try the development versions of
      \pkg{vegan} in
      \href{http://r-forge.r-project.org/projects/vegan/}{R-Forge} or
      \href{https://github.com/jarioksa/vegan}{GitHub} and report the
      problems and user experience to us.  }  } % end general

  \subsection{BUG FIXES}{

    \itemize{

      \item \code{envfit} function assumed that all external variables
      were either numeric or factors, and failed if they were, say,
      character strings. Now only numeric variables are taken as
      continuous vectors, and all other variables (character strings,
      logical) are coerced to factors if possible. The function also
      should work with degenerate data, like only one level of a
      factor or a constant value of a continuous environmental
      variable. The ties were wrongly in assessing permutation
      \eqn{P}-values in \code{vectorfit}.

      \item \code{nestednodf} with quantitative data was not
      consistent with binary models, and the fill was wrongly
      calculated with quantitative data.
      
      \item \code{oecosimu} now correctly adapts displayed quantiles
      of simulated values to the \code{alternative} test direction.

      \item \code{renyiaccum} plotting failed if only one level of
      diversity \code{scale} was used.
    
    }
  } % bug fixes

  \subsection{NEW FEATURES}{ 
    \itemize{ 

      \item The Kempton and Taylor algorithm was found unreliable in
      \code{fisherfit} and \code{fisher.alpha}, and now the estimation
      of Fisher \eqn{\alpha}{alpha} is only based on the number of
      species and the number of individuals.  The estimation of
      standard errors and profile confidence intervals also had to be
      scrapped.

      \item \code{renyiaccum}, \code{specaccum} and
      \code{tsallisaccum} functions gained \code{subset} argument.

      \item \code{renyiaccum} can now add a \code{collector} curve to
      to the analysis. The collector curve is the diversity
      accumulation in the order of the sampling units. With an
      interesting ordering or sampling units this allows comparing
      actual species accumulations with the expected randomized
      accumulation.
      
      \item \code{specaccum} can now perform weighted accumulation
      using the sampling effort as weights.

      }

  } % new features 


} % end 2.0-10

\section{Changes in version 2.0-9}{

\itemize{

   \item This version is released due to changes in programming
   interface and testing procedures in \R{} 3.0.2. If you are using an
   older version of \R, there is no need to upgrade \pkg{vegan}. There
   are no new features nor bug fixes. The only user-visible changes
   are in documentation and in output messages and formatting. Because
   of \R{} changes, this version is dependent on \R{} version 2.14.0
   or newer and on \pkg{lattice} package.

   }
}

\section{Changes in version 2.0-8}{

  \subsection{GENERAL}{
     \itemize{

       \item This is a maintenance release that fixes some issues
       raised by changed in \R{} toolset for processing vignettes.  In
       the same we also fix some typographic issues in the vignettes.

     }
  } % general

  \subsection{NEW FEATURES}{
    \itemize{

      \item \code{ordisurf} gained new arguments for more flexible
      definition of fitted models to better utilize the
      \pkg{mgcv}\code{::gam} function. 

      The linewidth of contours can
      now be set with the argument \code{lwd}.

      \item Labels to arrows are positioned in a better way in
      \code{plot} functions for the results of \code{envfit},
      \code{cca}, \code{rda} and \code{capscale}. The labels should no
      longer overlap the arrow tips.

      \item The setting test direction is clearer in \code{oecosimu}.

      \item \code{ordipointlabel} gained a \code{plot} method that can
      be used to replot the saved result.

    }
  } % new features

}

\section{Changes in version 2.0-7}{

  \subsection{NEW FUNCTIONS}{
    \itemize{

      \item \code{tabasco()} is a new function for graphical display
      of community data matrix.  Technically it is an interface to \R
      \code{heatmap}, but its use is closer to \pkg{vegan} function
      \code{vegemite}. The function can reorder the community data
      matrix similarly as \code{vegemite}, for instance, by ordination
      results. Unlike \code{heatmap}, it only displays dendrograms if
      supplied by the user, and it defaults to re-order the
      dendrograms by correspondence analysis.  Species are ordered to
      match site ordering or like determined by the user.

    }
  } % new functions

  \subsection{BUG FIXES}{
    \itemize{

      \item Function \code{fitspecaccum(..., model = "asymp")} fitted
      logistic model instead of asymptotic model (or the same as 
      \code{model = "logis"}).

      \item \code{nestedtemp()} failed with very sparse data (fill
      \eqn{< 0.38}\%).

    }
  } % bug fixes 

  \subsection{NEW FEATURES}{
    \itemize{

      \item The \code{plot} function for constrained ordination
      results (\code{cca}, \code{rda}, \code{capscale}) gained
      argument \code{axis.bp} (defaults \code{TRUE}) which can be used
      to suppress axis scale for biplot arrays.

      \item Number of iterations in nonmetric multidimensional scaling
      (NMDS) can be set with keyword \code{maxit} (defaults
      \code{200}) in \code{metaMDS}.

    }
  } % new features

  \subsection{DEPRECATED}{

    \itemize{

      \item The result objects of \code{cca}, \code{rda} and
      \code{capscale} will no longer have scores \code{u.eig},
      \code{v.eig} and \code{wa.eig} in the future versions of
      \pkg{vegan}. This change does not influence normal usage,
      because \pkg{vegan} functions do not need these items.  However,
      external scripts and packages may need changes in the future
      versions of \pkg{vegan}.

    }   
  } % deprecated

} % vegan 2.0-7

\section{Changes in version 2.0-6}{

  \subsection{BUG FIXES}{

    \itemize{

      \item The species scores were scaled wrongly in
      \code{capscale()}. They were scaled correctly only when Euclidean
      distances were used, but usually \code{capscale()} is used with
      non-Euclidean distances.  Most graphics will change and should be
      redone.  The change of scaling mainly influences the spread of
      species scores with respect to the site scores.

      \item Function \code{clamtest()} failed to set the minimum
      abundance threshold in some cases. In addition, the output was
      wrong when some of the possible species groups were missing. Both
      problems were reported by Richard Telford (Bergen, Norway).

      \item Plotting an object fitted by \code{envfit()} would fail if
      \code{p.max} was used and there were unused levels for one or
      more factors. The unused levels could result from deletion of
      observations with missing values or simply as the result of
      supplying a subset of a larger data set to \code{envfit()}.

      \item \code{multipart()} printed wrong information about the
      analysis type (but did the analysis correctly). Reported by
      Valerie Coudrain.

      \item \code{oecosimu()} failed if its \code{nestedfun} returned a
      data frame. A more fundamental fix will be in \pkg{vegan} 2.2-0,
      where the structure of the \code{oecosimu()} result will change.

      \item The plot of two-dimensional \code{procrustes()} solutions
      often draw original axes in a wrong angle. The problem was
      reported by Elizabeth Ottesen (MIT).

      \item Function \code{treedive()} for functional or phylogenetic
      diversity did not correctly match the species names between the
      community data and species tree when the tree contained species
      that did not occur in the data.  Related function
      \code{treedist()} for phylogenetic distances did not try to match
      the names at all.

    }
 
  } % bug fixes

  \subsection{NEW FEATURES}{

    \itemize{

      \item The output of \code{capscale()} displays the value of the
      additive constant when argument \code{add = TRUE} was used.
      
      \item \code{fitted()} functions for \code{cca()}, \code{rda()} and
      \code{capscale()} can now return conditioned (partial) component
      of the response: Argument \code{model} gained a new alternative
      \code{model = "pCCA"}. 

      \item \code{dispindmorisita()} output gained a new column for
      Chi-squared based probabilities that the null hypothesis (random
      distribution) is true.

      \item \code{metaMDS()} and \code{monoMDS()} have new default
      convergence criteria. Most importantly, scale factor of the
      gradient (\code{sfgrmin}) is stricter. The former limit was too
      slack with large data sets and iterations stopped early without
      getting close to the solution. In addition, \code{scores()}
      ignore now requests to dimensions beyond those calculated
      instead of failing, and \code{scores()} for \code{metaMDS()}
      results do not drop dimensions.

      \item \code{msoplot()} gained \code{legend} argument for
      positioning the legend.

      \item Nestedness function \code{nestednodf()} gained a \code{plot}
      method.

      \item \code{ordiR2step()} gained new argument \code{R2scope}
      (defaults \code{TRUE}) which can be used to turn off the criterion
      of stopping when the adjusted \eqn{R^2}{R-squared} of the current
      model exceeds that of the scope.  This option allows model
      building when the \code{scope} would be overdetermined (number of
      predictors higher than number of observations).

      \code{ordiR2step()} now handles partial redundancy analysis
      (pRDA).

      \item \code{orditorp()} gained argument \code{select} to select
      the rows or columns of the results to display.

      \item \code{protest()} prints the standardized residual statistic
      \eqn{m_{12}^2}{squared m12} in addition to the squared Procrustes
      correlation \eqn{R^2}{R-squared}.  Both were calculated, but only
      the latter was displayed.

      Permutation tests are much faster in \code{protest()}. Instead
      of calling repeatedly \code{procrustes()}, the goodness of fit
      statistic is evaluated within the function. 

      \item \code{wcmdscale()} gained methods for \code{print},
      \code{plot} etc. of the results. These methods are only used if
      the full \code{wcmdscale} result is returned with, e.g., argument
      \code{eig = TRUE}. The default is still to return only a matrix of
      scores similarly as the standard \R function \code{cmdscale()},
      and in that case the new methods are not used.

    }

  } % new features
   

} % end 2.0-6

\section{Changes in version 2.0-5}{

  \subsection{BUG FIXES}{

    \itemize{

      \item \code{anova(<cca_object>, ...)} failed with
      \code{by = "axis"} and \code{by = "term"}. The bug was reported by
      Dr Sven Neulinger (Christian Albrecht University, Kiel, Germany).

      \item \code{radlattice} did not honour argument \code{BIC = TRUE},
      but always displayed AIC.

    }

  } % bug fixes

  \subsection{NEW FUNCTIONS}{
    \itemize{

      \item Most \pkg{vegan} functions with permutation tests have now a
      \code{density} method that can be used to find empirical
      probability distributions of permutations. There is a new
      \code{plot} method for these functions that displays both the
      density and the observed statistic. The \code{density} function is
      available for \code{adonis}, \code{anosim}, \code{mantel},
      \code{mantel.partial}, \code{mrpp}, \code{permutest.cca} and
      \code{procrustes}.

      Function \code{adonis} can return several statistics, and it has
      now a \code{densityplot} method (based on \pkg{lattice}).

      Function \code{oecosimu} already had \code{density} and
      \code{densityplot}, but they are now similar to other \pkg{vegan}
      methods, and also work with \code{adipart}, \code{hiersimu} and
      \code{multipart}.

      \item \code{radfit} functions got a \code{predict} method that
      also accepts arguments \code{newdata} and \code{total} for new
      ranks and site totals for prediction.  The functions can also
      interpolate to non-integer \dQuote{ranks}, and in some models
      also extrapolate.

    }
  } % new functions

  \subsection{NEW FEATURES}{
    \itemize{

      \item Labels can now be set in the \code{plot} of \code{envfit}
      results. The labels must be given in the same order that the
      function uses internally, and new support function \code{labels}
      can be used to display the default labels in their correct order.

      \item Mantel tests (functions \code{mantel} and
      \code{mantel.partial}) gained argument \code{na.rm} which can be
      used to remove missing values. This options should be used with
      care: Permutation tests can be biased if the missing values were
      originally in matching or fixed positions.

      \item \code{radfit} results can be consistently accessed with
      the same methods whether they were a single model for a single
      site, all models for a single site or all models for all sites
      in the data.  All functions now have methods \code{AIC},
      \code{coef}, \code{deviance}, \code{logLik}, \code{fitted},
      \code{predict} and \code{residuals}.

    }

  } % new features
  
  \subsection{INSTALLATION AND BUILDING}{
    \itemize{

      \item Building of \pkg{vegan} vignettes failed with the latest
      version of LaTeX (TeXLive 2012).

      \item \R{} versions later than 2.15-1 (including development
      version) report warnings and errors when installing and checking
      \pkg{vegan}, and you must upgrade \pkg{vegan} to this version.
      The warnings concern functions \code{cIndexKM} and
      \code{betadisper}, and the error occurs in \code{betadisper}.
      These errors and warnings were triggered by internal changes in
      \R.
      
    }
  } % installation and building
  

} % version 2.0-5


\section{Changes in version 2.0-4}{

  \subsection{BUG FIXES}{
    \itemize{
    
      \item \code{adipart} assumed constant gamma diversity in
        simulations when assessing the \eqn{P}-value.  This could give
        biased results if the null model produces variable gamma
        diversities and option \code{weights = "prop"} is used.  The
        default null model (\code{"r2dtable"}) and the default option
        (\code{weights = "unif"}) were analysed correctly.

      \item \code{anova(<prc-object>, by = "axis")} and other
        \code{by} cases failed due to \file{NAMESPACE} issues.

      \item \code{clamtest} wrongly used frequencies instead of the
        counts when calculating sample coverage.  No detectable
        differences were produced when rerunning examples from Chazdon
        et al. 2011 and \pkg{vegan} help page.

      \item \code{envfit} failed with unused factor levels.

      \item \code{predict} for \code{cca} results with
        \code{type = "response"} or \code{type = "working"} failed with
        \code{newdata} if the number of rows did not match with the
        original data.  Now the \code{newdata} is ignored if it has a
        wrong number of rows.  The number of rows must match because
        the results in \code{cca} must be weighted by original row
        totals.  The problem did not concern \code{rda} or
        \code{capscale} results which do not need row weights.
        Reported by Glenn De'ath.
    
    }
  }% end bug fixes

  \subsection{NEW FEATURES}{
    \itemize{

      \item Functions for diversity partitioning (\code{adipart},
        \code{hiersimu} and \code{multipart}) have now \code{formula}
        and \code{default} methods.  The \code{formula} method is
        identical to the previous functions, but the \code{default}
        method can take two matrices as input. 

 	Functions \code{adipart} and \code{multipart} can be used for
 	fast and easy overall partitioning to alpha, beta and gamma
 	diversities by omitting the argument describing the hierarchy.

      \item The method in \code{betadisper} is biased with small
        sample sizes.  The effects of the bias are strongest with
        unequal sample sizes. A bias adjusted version was developed by
        Adrian Stier and Ben Bolker, and can be invoked with argument
        \code{bias.adjust} (defaults to \code{FALSE}).

      \item \code{bioenv} accepts dissimilarities (or square matrices
        that can be interpreted as dissimilarities) as an alternative to
        community data.  This allows using other dissimilarities than
        those available in \code{vegdist}.

      \item \code{plot} function for \code{envfit} results gained new
        argument \code{bg} that can be used to set background colour for
        plotted labels.

      \item \code{msoplot} is more configurable, and allows, for
        instance, setting y-axis limits.

      \item Hulls and ellipses are now filled using semitransparent
        colours in \code{ordihull} and \code{ordiellipse}, and the
        user can set the degree of transparency with a new argument
        \code{alpha}.  The filled shapes are used when these functions
        are called with argument \code{draw = "polygon"}.  Function
        \code{ordihull} puts labels (with argument \code{label = TRUE})
        now in the real polygon centre.

      \item \code{ordiplot3d} returns function \code{envfit.convert}
        and the projected location of the \code{origin}.  Together
        these can be used to add \code{envfit} results to existing
        \code{ordiplot3d} plots.

        Equal aspect ratio cannot be set exactly in \code{ordiplot3d}
        because underlying core routines do not allow this. Now
        \code{ordiplot3d} sets equal axis ranges, and the documents
        urge users to verify that the aspect ratio is reasonably equal
        and the graph looks like a cube.  If the problems cannot be
        solved in the future, \code{ordiplot3d} may be removed from
        next releases of \pkg{vegan}.

      \item Function \code{ordipointlabel} gained argument to
        \code{select} only some of the items for plotting. The
        argument can be used only with one set of points.

    }
  } % end new features 	 

}%end version 2.0-4

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
      by Eduard Szöcs (Uni Landau, Germany).

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
     step-across (extended) dissimilarities as
     \code{metaMDS(..., engine = "isoMDS")}. The change of defaults can
     also influence triggering of step-across in
     \code{capscale(..., metaMDSdist = TRUE)}.

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
   	  

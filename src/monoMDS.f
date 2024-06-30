      SUBROUTINE monoMDS (NOBJ, NFIX, NDIM, NDIS, NGRP,
     .  DISS, IIDX, JIDX, XINIT, ISTART,
     .  ISFORM, ITIES, IREGN, ISCAL, MAXITS, SRATMX,
     .  STRMIN, SFGRMN, 
     .  DIST, DHAT, X, STRESS, STRS, ITERS, ICAUSE)
C
C Subroutine for multidimensional scaling.
C
C 1.00 March 28, 2011
C 1.01 April 6, 2011  - added argument STRS(NGRP) to return the stress
C                       for each of the NGRP groups of dissimilarities
C                       i.e., from each separate regression.
C 1.02 January 7, 2016 - fixed bug in MONREG so that huge 
C                        dissimilarities are correctly handled in
C                        creating the initial partition for monotone 
C                        regression with primary tie treatment.
C 1.03 April 11, 2023 - changed calculation of M in ASORT4 from
C                       exponentiation to use of ISHIFT to avoid a
C                       problem with some optimizing compilers
C 1.04 March 8, 2024  - removed the use of ALLOCATABLE arrays by making
C                       IWORK, GRAD, and GRLAST permananet local arrays
C
C Written by Dr. Peter R. Minchin
C            Department of Biological Sciences
C            Southern Illinois University Edwardsville
C            PO Box 1651
C            Edwardsville, IL 62026-1651, U.S.A.
C            Phone: +1-618-650-2975   FAX: +1-618-650-3174
C            Email: pminchi@siue.edu
C
C Starting from a supplied initial configuarion, uses steepest descent
C   to minimize Kruskal's stress, a measure of badness-of-fit of one
C   or more regressions of distances onto the supplied dissimilarities.
C
C This routine is very flexible in the form of the dissimilarity
C   matrices it can handle because dissimilarities are input as a 
C   vector (DISS) with the associated i- and j-indices in parallel
C   vectors (IIDX, JIDX).  There is no need for missing values, as
C   these index vectors specify which pair of points each dissimilarity
C   belongs to and missing dissimilarities are simply omitted.  Another
C   element of flexibility is that there can be more than one group of
C   dissimilarities and a separate regression of distance on 
C   dissimilarity is done for each group.  Local or "row conditional"
C   scaling can be performed by inputting each row of the dissimilarity
C   matrix (minus the diagonals) as a separate group.  Hybrid scaling
C   of Faith et. al (1987) can be done by inputting the complete 
C   dissimilarity matrix as group 1 and only those dissimilarities
C   below a specified threshold as group 2 (and also setting input
C   parameter IREGN to 3).  New points can be inserted into an
C   existing ordination by inputting dissimilarities from the
C   rectangular matrix in which the existing points are rows and the
C   new points are columns.  The parameter NFIX is then set to the
C   number of existing points (rows) and these will be the first 
C   NFIX points in XINIT, which supplies their coordinates.  Normally,
C   when NFIX>0, the scaling parameter, ISCAL, should be 0, so that 
C   the coordinates of the fixed points remain unchanged in the 
C   output ordination scores, X.
C
C========INPUT ARGUMENTS:
C
C NOBJ = number of objects being ordinated
C
C NFIX = number of fixed objects (< NOBJ - if NFIX>0 the input
C        dissim matrix is assumed to be rectangular, with the
C        first NFIX objects, rows, remaining fixed and the rest
C        able to move)
C
C NDIM = number of dimensions for the MDS ordination
C
C NDIS = total number of dissimilarities among objects
C
C NGRP = number of groups of dissimilarities (a separate regression
C        of distances on dissimilarities is performed for each
C        group)
C
C DISS(NDIS) = vector of dissimilarities among objects
C
C IIDX(NDIS) = i-indices of objects in the dissimilarity vector
C
C JIDX(NDIS) = j-indices of objects in the dissimilarity vector
C
C XINIT(NOBJ,NDIM) = initial coordinates for NOBJ objects on NDIM
C                    dimensions
C                     
C ISTART(NGRP) = index of starting position of each group of
C                dissimilarities within vectors DISS, IIDX and JIDX
C
C ISFORM = Kruskal stress formula to be used: 1 or 2
C
C ITIES = Kruskal tie treatment: 1 = weak/primary, 2 = strong/secondary
C
C IREGN = Regression type: 1 = Kruskal monotone, 2 = linear, 3 = hybrid
C         (Kruskal monotone for the first NGRP/2 groups, linear for
C         the remainder of the groups - when this is selected it is
C         assumed that NGRP>1 and NGRP is even)
C
C ISCAL = scaling of final ordination: 0 = preserved (output raw scores
C         with no standardization), 1 = normalize (center and adjust
C         RMS distance of points to centroid to 1.0)
C
C MAXITS = maximum number of iterations - at least 200 recommended
C
C SRATMX = maximum average stress ratio (iterations stop when average
C          stress ratio exceeds this value but is still < 1) - a value
C          of 0.9999 or higher is recommended
C
C STRMIN = minimum stress to stop (iterations stop if stress falls
C          below this value) - a value of 0.001 or lower is
C          recommended
C
C SFGRMN = minimum scale factor of the gradient (iterations stop if
C          the scale factor of the gradient drops below this)
C
C========OUTPUT ARGUMENTS:
C
C DIST(NDIS) = distances among objects in the final ordination
C
C DHAT(NDIS) = fitted distances in final regression(s) of distance
C              on dissimilarity
C
C X(NOBJ,NDIM) = final ordination coordinates of NOBJ objects in NDIM
C                dimensions
C
C STRESS = final value of stress
C
C STRS(NGRP) = vector of stress values from each separate regression
C
C ITERS = number of iterations performed
C
C ICAUSE = reason for termination of iterations: 1 = max iterations
C          used, 2 = stress fell below STRMIN, 3 = stress ratio
C          exceeded SRATMX, 4 = scale factor of gradient fell below
C          SFGRMN
C
C---INPUT ARGUMENTS
      INTEGER, INTENT(IN) :: NOBJ, NFIX, NDIM, NDIS, NGRP,
     .  ISFORM, ITIES, IREGN, ISCAL, MAXITS
      INTEGER, INTENT(IN) :: IIDX(NDIS), JIDX(NDIS), ISTART(NGRP)
      DOUBLE PRECISION, INTENT(IN) :: XINIT(NOBJ,NDIM), DISS(NDIS),
     .  SRATMX, STRMIN, SFGRMN
C---OUTPUT ARGUMENTS
      INTEGER, INTENT(OUT) :: ITERS, ICAUSE
      DOUBLE PRECISION, INTENT(OUT) :: X(NOBJ,NDIM), DIST(NDIS), 
     .  DHAT(NDIS), STRESS, STRS(NGRP)
C-Removed the use of allocatable arrays (v 1.04)
C---TEMPORARY ARRAYS
      INTEGER :: IWORK(NDIS)
      DOUBLE PRECISION :: GRAD(NOBJ,NDIM), GRLAST(NOBJ,NDIM)
C
      DOUBLE PRECISION :: STRLST, SQRTN, SRATF1, SRATF2, FNGRP,
     .  STRINC, COSAV, ACOSAV, SRATAV, STEP, FNDIM, SFGR, SRATIO,
     .  SFACT, TFACT, DMEAN, CAGRGL, SFGLST, STPNEW, SSFACB, SSFACT,
     .  PARAM(2)
C-Removed the use of allocatable arrays (v 1.04)
C ALLOCATE THE TEMPORARY ARRAYS NEEDED
C 
C      ALLOCATE (IWORK(NDIS), GRAD(NOBJ,NDIM), GRLAST(NOBJ,NDIM))
C
C INITIALIZE SOME PARAMETERS
C
      SQRTN=SQRT(REAL(NDIS))
      SRATF1=0.5*(1.0+SRATMX)
      SRATF2=1.0-SRATF1
      FNGRP=REAL(NGRP)
      STRINC=1.1
      STRESS=1.0D0
      COSAV=0.0
      ACOSAV=0.0
      SRATAV=0.8
      STEP=0.0
      FNDIM=REAL(NDIM)
      SFGR=SQRTN
C
C PREPARE DISSIMILARITIES FOR USE IN REGRESSION(S)
C
      DO IGRP=1,NGRP
        I1=ISTART(IGRP)
        IF (IGRP.LT.NGRP) THEN
          N=ISTART(IGRP+1)-I1
        ELSE
          N=NDIS+1-I1
        ENDIF
        IF (N.GT.0) THEN
C
C PRE-SORT IF NECESSARY. NOT REQUIRED FOR LINEAR REGRESSION (IREGN=2)
C     NOR FOR DISSIM GROUPS > NGRP/2 IN HYBRID REGRESSION (IREGN=3)
C
          IF (IREGN.EQ.1.OR.(IREGN.EQ.3.AND.IGRP.LE.NGRP/2)) THEN
            CALL ASORT4 (DISS(I1),N,IIDX(I1),JIDX(I1))
          ENDIF
        ENDIF
      ENDDO
C
C COPY INITIAL CONFIGURATION TO CURRENT CONFIGURATION
C
      CALL MACOPY (XINIT,NOBJ,NOBJ,NDIM,X,NOBJ)
C
C INITALIZE GRADIENT (WILL BE USED AS THE FIRST "LAST GRADIENT")
C
      CALL MAINIT (GRAD,NOBJ,NDIM,NOBJ,SQRT(1.0/FNDIM))
C=======================================================================
C
C START OF ITERATION LOOP
C
      DO IT=0,MAXITS
C---Reset backup counter
        NBACK=0
C
C SAVE LAST STRESS VALUE AND SET STRESS BACK TO ZERO
C
C---Jump back to here when a backup if required
   10   STRLST=STRESS
        STRESS=0.0D0
C
C Normalize the current configuration (unless ISCAL=0)
C
        IF (ISCAL.GT.0) THEN
          CALL NRMCON (X,NOBJ,NDIM,NOBJ,SSFACT)
        ELSE
          SSFACT=1.0D0
        ENDIF
C
C MOVE CURRENT GRADIENT TO LAST GRADIENT (WITH APPROPRIATE SCALING)
C   AND SAVE LAST GRADIENT SCALING FACTOR
C
        CALL MACOPY (GRAD,NOBJ,NOBJ,NDIM,GRLAST,NOBJ)
        CALL MAMAS (GRLAST,NOBJ,NOBJ,NDIM,SSFACT)
        SFGLST=SFGR*SSFACT
C
C CLEAR CURRENT GRADIENT
C
        CALL MAINIT (GRAD,NOBJ,NDIM,NOBJ,0.0D0)
C
C COMPUTE DISTANCES
C
        CALL CLCDIS (X,NOBJ,NDIM,DIST,IIDX,JIDX,NDIS)
C-----------------------------------------------------------------------
C
C LOOP OVER THE NGRP GROUPS OF DISSIMILARITIES
C
        DO IGRP=1,NGRP
          I1=ISTART(IGRP)
          IF (IGRP.LT.NGRP) THEN
            N=ISTART(IGRP+1)-I1
          ELSE
            N=NDIS+1-I1
          ENDIF
          IF (N.GT.0) THEN
C
C PERFORM REGRESSION OF DISTANCES ON DISSIMILARITIES
C
            IF (IREGN.EQ.1.OR.(IREGN.EQ.3.AND.IGRP.LE.NGRP/2)) THEN
              CALL MONREG (DISS(I1),DIST(I1),DHAT(I1),IIDX(I1),
     .          JIDX(I1),IWORK(I1),N,ITIES)
            ELSEIF (IREGN.EQ.2.OR.(IREGN.EQ.3.AND.IGRP.GT.NGRP/2)) THEN
              CALL LINREG (DISS(I1),DIST(I1),DHAT(I1),N,PARAM)
            ENDIF
C
C COMPUTE AND ACCUMULATE STRESS
C
            CALL CLCSTR (DIST(I1),DHAT(I1),N,SFACT,TFACT,STRS(IGRP),
     .        ISFORM,DMEAN)
            STRESS=STRESS+SFACT/TFACT
C
C ACCUMULATE THE NEGATIVE GRADIENT
C
            CALL CLCGRD (X,GRAD,NOBJ,NDIM,DIST(I1),DHAT(I1),
     .        IIDX(I1),JIDX(I1),N,STRS(IGRP),SFACT,TFACT,ISFORM,DMEAN)
          ENDIF
        ENDDO
C
C END OF DISSIMILARITY GROUP LOOP
C
C-----------------------------------------------------------------------
C
C COMPUTE OVERALL STRESS AND ITS REDUCTION RATIO
C
        STRESS=SQRT(STRESS/FNGRP)
        IF (IT.EQ.0) THEN
          SRATIO=0.8D0
        ELSE
          SRATIO=STRESS/STRLST
        ENDIF
C
C TO PREVENT FIXED POINTS BEING MOVED, ZERO THEIR GRADIENT ELEMENTS
C
        IF (NFIX.GT.0) THEN
          DO I=1,NFIX
            DO IDIM=1,NDIM
              GRAD(I,IDIM)=0.0D0
            ENDDO
          ENDDO
        ENDIF
C
C CALCULATE SCALE FACTOR OF NEW GRADIENT AND ITS DIRECTION COSINE
C   WITH THE LAST GRADIENT
C
        CALL CLCSFA (GRAD,GRLAST,NOBJ,NDIM,NOBJ,SFGR,CAGRGL,
     .    SFGLST)
C
C IF STRESS HAS INCREASED APPRECIABLY, OR THE ANGLE BETWEEN THE NEW
C   GRADIENT AND THE LAST GRADIENT IS APPROACHING 180 DEGREES, REDUCE
C   THE STEP SIZE AND BACK UP ALONG THE LAST GRADIENT
C   (BACKUP PROCEDURE IS PERFORMED A MAXIMUM OF 3 TIMES PER ITERATION)
C
        IF ((SRATIO.GT.STRINC.OR.CAGRGL.LT.(-0.95)).AND.NBACK.LT.3)
     .    THEN
          STPNEW=STEP*0.1D0
          CALL BACKUP (X,GRAD,GRLAST,NOBJ,NDIM,NOBJ,NBACK,SSFACT,
     .      SSFACB,STRESS,STRLST,SFGR,SFGLST,STEP,STPNEW)
          GO TO 10
        ENDIF
C
C UPDATE AVERAGES OF STRESS RATIO AND PREVIOUS DIRECTION COSINES
C
        SRATAV=SRATIO**0.33334D0 * SRATAV**0.66666D0
        COSAV=CAGRGL*0.66666D0 + COSAV*0.33334D0
        ACOSAV=ABS(CAGRGL)*0.66666D0 + ACOSAV*0.33334D0
C
C COMPUTE NEW STEP SIZE
C
        CALL CLCSTP (STEP,IT,SFGR,STRESS,COSAV,ACOSAV,SRATIO,
     .    SRATAV)
C
C CHECK WHETHER OR NOT TO KEEP ITERATING
C
C IF STRESS IN SMALL ENOUGH, TERMINATE
C
        ITERS=IT
        IF (STRESS.LT.STRMIN) THEN
          ICAUSE=2
          EXIT
C
C IF STRESS IS DECREASING SLOWLY, TERMINATE.
C
        ELSEIF (ABS(SRATAV-SRATF1).LE.SRATF2.AND.
     .          ABS(SRATIO-SRATF1).LE.SRATF2) THEN
          ICAUSE=3
          EXIT
C
C IF SCALE FACTOR OF GRADIENT IS SMALL ENOUGH, TERMINATE
C
        ELSEIF (SFGR.LE.SFGRMN) THEN
          ICAUSE=4
          EXIT
C
C IF THIS IS THE FINAL ITERATION, TERMINATE
C
        ELSEIF (IT.EQ.MAXITS) THEN
          ICAUSE=1
          EXIT
        ENDIF
C
C COMPUTE NEW CONFIGURATION
C
        CALL NEWCON (X,GRAD,NOBJ,NDIM,NOBJ,STEP,SFGR)
      ENDDO
C
C END OF ITERATION LOOP
C
C=======================================================================
C
C-Removed the use of allocatable arrays (v 1.04)
C DEALLOCATE THE TEMPORARY ARRAYS AND RETURN
C
C      DEALLOCATE (IWORK, GRAD, GRLAST)
      RETURN
      END SUBROUTINE monoMDS

      SUBROUTINE ASORT4 (X,N,I1,I2)
C
C SORTS THE REAL*8 VECTOR X ASCENDING, SIMULTANEOUSLY PERMUTING THE
C   INTEGER VECTORS I1 AND I2 IN PARALLEL.
C
C ADAPTED FROM SUPERCHARGED SHELL SORT ALGORITHM PUBLISHED ON PAGE 488
C   OF THE MAY 1983 EDITION OF BYTE MAGAZINE.
C
C Adapted by Dr. Peter R. Minchin
C            Department of Biological Sciences
C            Southern Illinois University Edwardsville
C            PO Box 1651
C            Edwardsville, IL 62026-1651, U.S.A.
C            Phone: +1-618-650-2975   FAX: +1-618-650-3174
C            Email: pminchi@siue.edu
C
      INTEGER I1(N), I2(N), I1TEMP, I2TEMP
      DOUBLE PRECISION X(N), XTEMP
C
      IF (N.LT.2) RETURN
      FN=REAL(N)
      NLOOPS=MAX(NINT(LOG(FN)/LOG(2.)),1)
C REPLACED EXPONENTIATION BY ISHFT TO AVOID A PROBLEM WITH SOME
C   OPTIMIZING COMPILERS (v 1.03)
C      M=2**(NLOOPS-1)
      M=ISHFT(1, NLOOPS-1)
      DO II=1,NLOOPS
        FM=M
        DO I=1,MAX(1,N-M)
          IF (X(I).GT.X(I+M)) THEN
            XTEMP=X(I+M)
            I1TEMP=I1(I+M)
            I2TEMP=I2(I+M)
            X(I+M)=X(I)
            I1(I+M)=I1(I)
            I2(I+M)=I2(I)
            IF (I.LE.M) THEN
              X(I)=XTEMP
              I1(I)=I1TEMP
              I2(I)=I2TEMP
            ELSE
              DO J=I-M,1,-M
                IF (XTEMP.GE.X(J)) EXIT
                  X(J+M)=X(J)
                  I1(J+M)=I1(J)
                  I2(J+M)=I2(J)
              ENDDO
              X(J+M)=XTEMP
              I1(J+M)=I1TEMP
              I2(J+M)=I2TEMP
            ENDIF
          ENDIF
        ENDDO
        M=INT(FM/2.)
      ENDDO
      RETURN
      END SUBROUTINE ASORT4

      SUBROUTINE BACKUP (X,GRAD,GRLAST,NOBJ,NDIM,MAXOBJ,NBACK,SSFACT,
     .  SSFACB,STRESS,STRLST,SFGR,SFGLST,STEP,STPNEW)
C
C BACKS THE CONFIGURATION UP ALONG THE LAST GRADIENT
C
C Written by Dr. Peter R. Minchin
C            Department of Biological Sciences
C            Southern Illinois University Edwardsville
C            PO Box 1651
C            Edwardsville, IL 62026-1651, U.S.A.
C            Phone: +1-618-650-2975   FAX: +1-618-650-3174
C            Email: pminchi@siue.edu
C
      DOUBLE PRECISION X(MAXOBJ,NDIM), GRAD(MAXOBJ,NDIM), 
     .  GRLAST(MAXOBJ,NDIM), SSFACT, SSFACB, STRESS, STRLST, SFGR,
     .  SFGLST,STEP,STPNEW,FACTOR
C
      NBACK=NBACK+1
      IF (NBACK.EQ.1) THEN
        SSFACB=1.0D0
      ELSE
        SSFACB=SSFACB*SSFACT
      ENDIF
      FACTOR=SSFACB*(STEP-STPNEW)/SFGLST
      DO IDIM=1,NDIM
        DO IOBJ=1,NOBJ
          X(IOBJ,IDIM)=X(IOBJ,IDIM)-FACTOR*GRLAST(IOBJ,IDIM)
          GRAD(IOBJ,IDIM)=GRLAST(IOBJ,IDIM)
        ENDDO
      ENDDO
      STEP=STPNEW
      SFGR=SFGLST
      STRESS=STRLST
      RETURN
      END SUBROUTINE BACKUP

      SUBROUTINE CLCDIS (X,MAXOBJ,NDIM,DIST,IIDX,JIDX,NDIS)
C
C COMPUTES EUCLIDEAN DISTANCES BETWEEN EACH PAIR OF THE NOBJ POINTS
C   WHOSE CO-ORDINATES IN NDIM DIMENSIONS ARE IN X(NOBJ,NDIM).
C
C THE DISTANCES ARE STORED IN THE VECTOR DIST AND THE INDEX VECTORS
C   IIDX AND JIDX HOLD THE I,J INDICES OF THE POINTS IN THE ORDER IN
C   WHICH THE COMPARISONS ARE TO BE MADE.
C
C Written by Dr. Peter R. Minchin
C            Department of Biological Sciences
C            Southern Illinois University Edwardsville
C            PO Box 1651
C            Edwardsville, IL 62026-1651, U.S.A.
C            Phone: +1-618-650-2975   FAX: +1-618-650-3174
C            Email: pminchi@siue.edu
C
      INTEGER IIDX(NDIS), JIDX(NDIS)
      DOUBLE PRECISION X(MAXOBJ,NDIM), DIST(NDIS)
C
C INITIALIZE DISTANCES
C
      DO I=1,NDIS
        DIST(I)=0.0
      ENDDO
C
C ACCUMULATE SUMS OF SQUARED DIFFERENCES ON EACH DIMENSION
C
      DO IDIM=1,NDIM
        DO K=1,NDIS
          DIST(K)=DIST(K)+(X(IIDX(K),IDIM)-X(JIDX(K),IDIM))**2
        ENDDO
      ENDDO
C
C TAKE SQUARE ROOTS OF TOTALS
C
      DO I=1,NDIS
        DIST(I)=SQRT(DIST(I))
      ENDDO
      RETURN
      END SUBROUTINE CLCDIS

      SUBROUTINE CLCGRD (X,GRAD,MAXOBJ,NDIM,DIST,DHAT,IIDX,JIDX,
     .  NDIS,STRESS,SFACT,TFACT,ISFORM,DMEAN)
C
C ACCUMULATES THE NEGATIVE GRADIENT IN GRAD(NOBJ,NDIM).
C
C Written by Dr. Peter R. Minchin
C            Department of Biological Sciences
C            Southern Illinois University Edwardsville
C            PO Box 1651
C            Edwardsville, IL 62026-1651, U.S.A.
C            Phone: +1-618-650-2975   FAX: +1-618-650-3174
C            Email: pminchi@siue.edu
C
      INTEGER IIDX(NDIS), JIDX(NDIS)
      DOUBLE PRECISION X(MAXOBJ,NDIM), GRAD(MAXOBJ,NDIM), DIST(NDIS),
     .  DHAT(NDIS), STRESS, SFACT, TFACT, DMEAN, SOTSQ, RECIPT,
     .  DELTA
C
      IF (STRESS.LE.0.0D0) RETURN
      SOTSQ=SFACT/(TFACT**2)
      RECIPT=1.0/TFACT
      DO IDIM=1,NDIM
        IF (ISFORM.LE.1) THEN
C---Kruskal's stress formula 1
          DO K=1,NDIS
            IF (DIST(K).GT.0.0) THEN
              DELTA=(SOTSQ-RECIPT*(DIST(K)-DHAT(K))/DIST(K))*
     .          (X(IIDX(K),IDIM)-X(JIDX(K),IDIM))
              GRAD(IIDX(K),IDIM)=GRAD(IIDX(K),IDIM)+DELTA
              GRAD(JIDX(K),IDIM)=GRAD(JIDX(K),IDIM)-DELTA
            ENDIF
          ENDDO
        ELSE
C---Kruskal's stress formula 2
          DO K=1,NDIS
            IF (DIST(K).GT.0.0) THEN
              DELTA=(SOTSQ*(DIST(K)-DMEAN)/DIST(K)-
     .          RECIPT*(DIST(K)-DHAT(K))/DIST(K))*
     .          (X(IIDX(K),IDIM)-X(JIDX(K),IDIM))
              GRAD(IIDX(K),IDIM)=GRAD(IIDX(K),IDIM)+DELTA
              GRAD(JIDX(K),IDIM)=GRAD(JIDX(K),IDIM)-DELTA
            ENDIF
          ENDDO
        ENDIF
      ENDDO
      RETURN
      END SUBROUTINE CLCGRD

      SUBROUTINE CLCSFA (GRAD,GRLAST,NOBJ,NDIM,MAXOBJ,SFGR,CAGRGL,
     .  SFGLST)
C
C COMPUTES SCALE FACTOR OF THE NEGATIVE GRADIENT IN GRAD(NOBJ,NDIM)
C   AND ALSO FINDS ITS DIRECTION COSINE WITH THE PREVIOUS GRADIENT,
C   STORED IN GRLAST(NOBJ,NDIM).
C
C Written by Dr. Peter R. Minchin
C            Department of Biological Sciences
C            Southern Illinois University Edwardsville
C            PO Box 1651
C            Edwardsville, IL 62026-1651, U.S.A.
C            Phone: +1-618-650-2975   FAX: +1-618-650-3174
C            Email: pminchi@siue.edu
C
      DOUBLE PRECISION GRAD(MAXOBJ,NDIM), GRLAST(MAXOBJ,NDIM),
     .  SFGR, CAGRGL, SFGLST, FN, FACTOR
C
      FN=DBLE(NOBJ)
      SFGR=0.0D0
      CAGRGL=0.0D0
      DO IDIM=1,NDIM
        DO IOBJ=1,NOBJ
          SFGR=SFGR+GRAD(IOBJ,IDIM)**2
          CAGRGL=CAGRGL+GRAD(IOBJ,IDIM)*GRLAST(IOBJ,IDIM)
        ENDDO
      ENDDO
      SFGR=SQRT(SFGR/FN)
      FACTOR=SFGR*SFGLST*FN
      IF (FACTOR.NE.0.0) CAGRGL=CAGRGL/FACTOR
      RETURN
      END SUBROUTINE CLCSFA

      SUBROUTINE CLCSTP (STEP,IT,SFGR,STRESS,COSAV,ACOSAV,SRATIO,
     .  SRATAV)
C
C COMPUTES NEW STEP SIZE.  BASED ON J.B. KRUSKAL'S STEP SIZE ALGORITHM
C   IN THE BELL LABORATORIES PROGRAM "KYST".
C
C Written by Dr. Peter R. Minchin
C            Department of Biological Sciences
C            Southern Illinois University Edwardsville
C            PO Box 1651
C            Edwardsville, IL 62026-1651, U.S.A.
C            Phone: +1-618-650-2975   FAX: +1-618-650-3174
C            Email: pminchi@siue.edu
C
      DOUBLE PRECISION STEP, SFGR, STRESS, COSAV, ACOSAV, SRATIO,
     .  SRATAV, FACTR1, FACTR2, FACTR3
C
C MAKE A WILD GUESS FOR INITIAL STEP SIZE (IT=0).
C
      IF (IT.EQ.0) THEN
        STEP=25.0*STRESS*SFGR
C
C OTHERWISE, UPDATE STEP SIZE.
C
      ELSE
        FACTR1=4.0**COSAV
        FACTR2=1.6/( (1.0+(MIN(1D0,SRATAV))**5) *
     .    (1.0+ACOSAV-ABS(COSAV)) )
        FACTR3=SQRT(MIN(1D0,SRATIO))
        STEP=STEP*FACTR1*FACTR2*FACTR3
      ENDIF
      RETURN
      END SUBROUTINE CLCSTP

      SUBROUTINE CLCSTR (DIST,DHAT,N,SNUMER,SDENOM,STRESS,ISFORM,DMEAN)
C
C COMPUTES KRUSKAL'S STRESS (FORMULA 1 or 2)
C
C Written by Dr. Peter R. Minchin
C            Department of Biological Sciences
C            Southern Illinois University Edwardsville
C            PO Box 1651
C            Edwardsville, IL 62026-1651, U.S.A.
C            Phone: +1-618-650-2975   FAX: +1-618-650-3174
C            Email: pminchi@siue.edu
C
      DOUBLE PRECISION DIST(N), DHAT(N), SNUMER, SDENOM, STRESS,
     .  DMEAN
C
      SNUMER=0.0D0
      SDENOM=0.0D0
      DMEAN=0.0D0
      IF (ISFORM.GE.2) THEN
        DO I=1,N
          DMEAN=DMEAN+DIST(I)
        ENDDO
        DMEAN=DMEAN/DBLE(N)
        DO I=1,N
          SNUMER=SNUMER+(DIST(I)-DHAT(I))**2
          SDENOM=SDENOM+(DIST(I)-DMEAN)**2
        ENDDO
      ELSE
        DO I=1,N
          SNUMER=SNUMER+(DIST(I)-DHAT(I))**2
          SDENOM=SDENOM+DIST(I)**2
        ENDDO
      ENDIF
      STRESS=SQRT(SNUMER/SDENOM)
      RETURN
      END SUBROUTINE CLCSTR

      SUBROUTINE LINREG (DISS,DIST,DHAT,N,C)
C
C PERFORMS LINEAR REGRESSION OF DIST ON DISS, PLACING THE FITTED VALUES
C   IN DHAT.  THE REGRESSION COEFFICIENTS ARE RETURNED IN C:
C   C(1) = INTERCEPT, C(2) = SLOPE.
C
C Written by Dr. Peter R. Minchin
C            Department of Biological Sciences
C            Southern Illinois University Edwardsville
C            PO Box 1651
C            Edwardsville, IL 62026-1651, U.S.A.
C            Phone: +1-618-650-2975   FAX: +1-618-650-3174
C            Email: pminchi@siue.edu
C
      DOUBLE PRECISION DISS(N), DIST(N), DHAT(N),
     .  AVDISS, SSDISS, AVDIST, SPROD, C(2), FN
C
C COMPUTE MEANS, SUMS OF SQUARES AND SUM OF PRODUCTS
C
      FN=DBLE(N)
      AVDIST=0.0D0
      AVDISS=0.0D0
      DO K=1,N
        AVDIST=AVDIST+DIST(K)
        AVDISS=AVDISS+DISS(K)
      ENDDO
      AVDIST=AVDIST/FN
      AVDISS=AVDISS/FN
      SSDISS=0.0D0
      SPROD=0.0D0
      DO K=1,N
        SSDISS=SSDISS+(DISS(K)-AVDISS)**2
        SPROD=SPROD+(DIST(K)-AVDIST)*(DISS(K)-AVDISS)
      ENDDO
C
C COMPUTE SLOPE AND INTERCEPT
C
      C(2)=SPROD/SSDISS
      C(1)=AVDIST-C(2)*AVDISS
C
C CALCULATE FITTED VALUES
C
      DO K=1,N
        DHAT(K)=C(1)+C(2)*DISS(K)
      ENDDO
      RETURN
      END SUBROUTINE LINREG

      SUBROUTINE MACOPY (A,MAXN1,N,M,B,MAXN2)
C
C COPY MATRIX A(N,M) TO B(N,M).
C
C CAN BE USED TO PACK AN ARRAY BY PASSING THE SAME ARRAY AS A AND B,
C   WITH MAXN2 SET TO N.
C
C Written by Dr. Peter R. Minchin
C            Department of Biological Sciences
C            Southern Illinois University Edwardsville
C            PO Box 1651
C            Edwardsville, IL 62026-1651, U.S.A.
C            Phone: +1-618-650-2975   FAX: +1-618-650-3174
C            Email: pminchi@siue.edu
C
      DOUBLE PRECISION A(MAXN1,M), B(MAXN2,M)
C
      DO J=1,M
        DO I=1,N
          B(I,J)=A(I,J)
        ENDDO
      ENDDO
      RETURN
      END SUBROUTINE MACOPY

      SUBROUTINE MAINIT (X,M,N,MAXM,CONST)
C
C INITIALIZES MATRIX X(M,N) WITH THE VALUE CONST.
C
C Written by Dr. Peter R. Minchin
C            Department of Biological Sciences
C            Southern Illinois University Edwardsville
C            PO Box 1651
C            Edwardsville, IL 62026-1651, U.S.A.
C            Phone: +1-618-650-2975   FAX: +1-618-650-3174
C            Email: pminchi@siue.edu
C
      DOUBLE PRECISION X(MAXM,N), CONST
C
      DO J=1,N
        DO I=1,M
          X(I,J)=CONST
        ENDDO
      ENDDO
      RETURN
      END SUBROUTINE MAINIT

      SUBROUTINE MAMAS (A,MAXL,L,M,S)
C
C MULTIPLY MATRIX A(L,M) BY THE SCALAR S
C
C Written by Dr. Peter R. Minchin
C            Department of Biological Sciences
C            Southern Illinois University Edwardsville
C            PO Box 1651
C            Edwardsville, IL 62026-1651, U.S.A.
C            Phone: +1-618-650-2975   FAX: +1-618-650-3174
C            Email: pminchi@siue.edu
C
      DOUBLE PRECISION A(MAXL,M), S
C
      DO I=1,L
        DO J=1,M
          A(I,J)=A(I,J)*S
        ENDDO
      ENDDO
      RETURN
      END SUBROUTINE MAMAS

      SUBROUTINE MONREG (DISS,DIST,DHAT,IIDX,JIDX,IWORK,N,ITIES)
C
C PERFORMS KRUSKAL'S MONOTONE REGRESSION OF DIST ON DISS, PLACING THE 
C   FITTED VALUES IN DHAT.
C
C ASSUMES THAT DISS HAS BEEN SORTED ASCENDING.
C
C Written by Dr. Peter R. Minchin
C            Department of Biological Sciences
C            Southern Illinois University Edwardsville
C            PO Box 1651
C            Edwardsville, IL 62026-1651, U.S.A.
C            Phone: +1-618-650-2975   FAX: +1-618-650-3174
C            Email: pminchi@siue.edu
C
      INTEGER IIDX(N), JIDX(N)
      INTEGER IWORK(N)
      DOUBLE PRECISION DISS(N), DIST(N), DHAT(N), DNEXT, TOLER, SUM,
     .  DHATAV
C---Tolerance used to test if dissimilarities are equal
      TOLER=SQRT(EPSILON(SUM))
C***********************************************************************
C
C CREATE INITIAL PARTITION:
C
C   PRE-PROCESS BLOCKS OF TIED DISS VALUES DEPENDING ON THE SELECTED
C     TIE TREATMENT (ITIES)
C-----------------------------------------------------------------------
      J=0
      DO I=1,N
        IF (I.LT.N) THEN
          DNEXT=DISS(I+1)
        ELSE
C---Bug fix January 7, 2016: correctly handles huge dissimilarities
C          DNEXT=DISS(I)+1.0
          DNEXT=DISS(I)*2.0
        END IF
        IF (ABS(DNEXT-DISS(I)).GT.TOLER) THEN
C---NTIE is the number of DISS values in the current group of tied values
          NTIE=I-J
          IF (NTIE.GT.1) THEN
            IF (ITIES.LE.1) THEN
C---Primary tie treatment: sort DIST values within this tied group into
C     ascending order, permuting the index vectors accordingly.
C     Initialize fitted values, DHAT, to be equal to DIST and make each
C     value an intial block of size 1.
              CALL ASORT4 (DIST(J+1),NTIE,IIDX(J+1),JIDX(J+1))
              DO K=J+1,I
                DHAT(K)=DIST(K)
                IWORK(K)=1
              ENDDO
            ELSE
C---Secondary tie treatment: set fitted value for all members of this
C     tied group to the average of the DIST values for the group.  Make
C     the group of tied values an initial block equal in size to the
C     number of tied values.
              IF (NTIE.EQ.2) THEN
                SUM=DIST(J+1)+DIST(I)
              ELSE
                SUM=0.0
                DO K=J+1,I
                  SUM=SUM+DIST(K)
                ENDDO
              ENDIF
              DHAT(J+1)=SUM
              DHAT(I)=SUM
              IWORK(J+1)=NTIE
              IWORK(I)=NTIE
            ENDIF
          ELSE
C---Initially, make each unique (not tied) value is DISS a block of
C   size 1.
            IWORK(I)=1
            DHAT(I)=DIST(I)
          ENDIF
          J=I
        ENDIF
      ENDDO
C***********************************************************************
C
C NOW START THE MONOTONE REGRESSION PROCEDURE
C
C ICURR is the block we are currently examining.
C
C WE START BY LOOKING AT THE FIRST BLOCK.
C-----------------------------------------------------------------------
      ICURR=1
C-----------------------------------------------------------------------
C START OF PROCEDURE FOR THE CURRENT BLOCK.  IT IS INITIALLY UP-ACTIVE,
C   AND NEITHER UP- NOR DOWN-SATISIFIED.
C
C UP-ACTIVE  (IACTIV=1) means we are comparing DHAT for the current
C   block with DHAT for blocks to its right.
C
C DOWN-ACTIVE (IACTIV=0) means we are comparing DHAT for the current
C   block with DHAT for blocks to its left.
C
C UP-SATISFIED means that DHAT for the current block is less than DHAT
C   for the block to its right.
C
C DOWN-SATISFIED means that DHAT for the current block is greater than
C   DHAT for the block to its left.
C
C NSATIS will be 2 once a block is both UP- and DOWN-SATISFIED.
C
C-----------------------------------------------------------------------
   25 IACTIV=1
      NSATIS=0
C---Compute DHAT for the current block
      DHATAV=DHAT(ICURR)/IWORK(ICURR)
C-----------------------------------------------------------------------
C   ACCORDING TO CURRENT ACTIVITY, CHECK WHETHER CURRENT BLOCK IS
C     UP-SATISFIED (IACTIV=1) OR DOWN-SATISIFED (IACTIV=0)
C-----------------------------------------------------------------------
   30 IF (IACTIV.EQ.0) THEN
C-----------------------------------------------------------------------
C     CHECK WHETHER THIS BLOCK IS DOWN-SATISIFIED.  IF NOT, MERGE.
C-----------------------------------------------------------------------
        IF (ICURR.EQ.1) THEN
C---If it's the first block, it is by definition down-satisfied
          NSATIS=NSATIS+1
        ELSEIF (DHATAV.GT.DHAT(ICURR-1)/IWORK(ICURR-1)) THEN
C---Current block is down-satisfied
          NSATIS=NSATIS+1
        ELSE
C---Current block is not down-satisfied, so merge it with the block to
C   its left and make the new merged block the current block
          ICNEW=ICURR-IWORK(ICURR-1)
          JCNEW=ICURR+IWORK(ICURR)-1
          IWORK(ICNEW)=JCNEW-ICNEW+1
          IWORK(JCNEW)=IWORK(ICNEW)
          DHAT(ICNEW)=DHAT(ICURR-1)+DHAT(ICURR)
          DHAT(JCNEW)=DHAT(ICNEW)
          ICURR=ICNEW
          DHATAV=DHAT(ICURR)/IWORK(ICURR)
          NSATIS=0
        ENDIF
      ELSE
C-----------------------------------------------------------------------
C     CHECK WHETHER THIS BLOCK IS UP-SATISIFIED.  IF NOT, MERGE.
C-----------------------------------------------------------------------
C---Index of first member of the block to the right
        INEXT=ICURR+IWORK(ICURR)
        IF (INEXT.GT.N) THEN
C---If current block is last block, it is, by definition up-satisified
          NSATIS=NSATIS+1
        ELSEIF (DHATAV.LT.DHAT(INEXT)/IWORK(INEXT)) THEN
C---Current block is up-satisfied
          NSATIS=NSATIS+1
        ELSE
C---Current block is not up-satisfied, so merge it with the block to
C   its right and make the new merged block the current block
          IWORK(ICURR)=IWORK(ICURR)+IWORK(INEXT)
          JCURR=ICURR+IWORK(ICURR)-1
          IWORK(JCURR)=IWORK(ICURR)
          DHAT(ICURR)=DHAT(ICURR)+DHAT(INEXT)
          DHAT(JCURR)=DHAT(ICURR)
          DHATAV=DHAT(ICURR)/IWORK(ICURR)
          NSATIS=0
        ENDIF
      ENDIF
C-----------------------------------------------------------------------
C TOGGLE ACTIVITY FOR THE CURRENT BLOCK.  CHECK IF IT IS BOTH UP- AND
C   DOWN-SATISFIED.  IF NOT, GO BACK AND KEEP TRYING.
C-----------------------------------------------------------------------
      IACTIV=1-IACTIV
      IF (NSATIS.LT.2) GO TO 30
C-----------------------------------------------------------------------
C ADVANCE TO THE NEXT BLOCK, CHECK IF FINISHED
C-----------------------------------------------------------------------
      ICURR=ICURR+IWORK(ICURR)
      IF (ICURR.LE.N) GO TO 25
C***********************************************************************
C CHECKING AND MERGING OF BLOCKS IS COMPLETED.  NOW PUT APPROPRIATE
C   AVERAGE VALUES IN ALL ELEMENTS OF DHAT
C-----------------------------------------------------------------------
      ICURR=1
   40   IF (IWORK(ICURR).GT.2) THEN
          DHATAV=DHAT(ICURR)/IWORK(ICURR)
          DO J=ICURR,ICURR+IWORK(ICURR)-1
            DHAT(J)=DHATAV
          ENDDO
        ELSEIF (IWORK(ICURR).EQ.2) THEN
          DHAT(ICURR)=0.5*DHAT(ICURR)
          DHAT(ICURR+1)=DHAT(ICURR)
        ENDIF
        ICURR=ICURR+IWORK(ICURR)
      IF (ICURR.LT.N) GO TO 40
      RETURN
      END SUBROUTINE MONREG

      SUBROUTINE NEWCON (X,GRAD,NOBJ,NDIM,MAXOBJ,STEP,SFGR)
C
C COMPUTES NEW CONFIGURATION, CHANGING EACH CO-ORDINATE ACCORDING TO
C   THE NEGATIVE GRADIENT.
C
C Written by Dr. Peter R. Minchin
C            Department of Biological Sciences
C            Southern Illinois University Edwardsville
C            PO Box 1651
C            Edwardsville, IL 62026-1651, U.S.A.
C            Phone: +1-618-650-2975   FAX: +1-618-650-3174
C            Email: pminchi@siue.edu
C
      DOUBLE PRECISION X(MAXOBJ,NDIM), GRAD(MAXOBJ,NDIM), STEP, SFGR,
     .  FACTOR
C
      FACTOR=STEP/SFGR
      DO IDIM=1,NDIM
        DO IOBJ=1,NOBJ
          X(IOBJ,IDIM)=X(IOBJ,IDIM)+FACTOR*GRAD(IOBJ,IDIM)
        ENDDO
      ENDDO
      RETURN
      END SUBROUTINE NEWCON

      SUBROUTINE NRMCON (X,NOBJ,NDIM,MAXOBJ,SSFACT)
C
C NORMALIZES THE CURRENT CONFIGURATION IN X(NOBJ,NDIM).  THE CENTROID
C   IS MOVED TO THE ZERO ORIGIN AND THE SUM OF SQUARES ADJUSTED SUCH
C   THAT THE RMS DISTANCE BETWEEN POINTS AND THE ORIGIN IS 1.0.
C
C ON RETURN, SSFACT CONTAINS THE SCALING FACTOR WHICH WAS USED.
C
C Written by Dr. Peter R. Minchin
C            Department of Biological Sciences
C            Southern Illinois University Edwardsville
C            PO Box 1651
C            Edwardsville, IL 62026-1651, U.S.A.
C            Phone: +1-618-650-2975   FAX: +1-618-650-3174
C            Email: pminchi@siue.edu
C
      DOUBLE PRECISION X(MAXOBJ,NDIM), SSFACT, FN, FMEAN
C
C CENTER THE CONFIGURATION AND COMPUTE SCALING FACTOR
C
      FN=DBLE(NOBJ)
      SSFACT=0.0D0
      DO IDIM=1,NDIM
        FMEAN=0.0D0
        DO IOBJ=1,NOBJ
          FMEAN=FMEAN+X(IOBJ,IDIM)
        ENDDO
        FMEAN=FMEAN/FN
        DO IOBJ=1,NOBJ
          X(IOBJ,IDIM)=X(IOBJ,IDIM)-FMEAN
          SSFACT=SSFACT+X(IOBJ,IDIM)**2
        ENDDO
      ENDDO
      SSFACT=SQRT(FN/SSFACT)
C
C SCALE THE CONFIGURATION
C
      CALL MAMAS (X,MAXOBJ,NOBJ,NDIM,SSFACT)
      RETURN
      END SUBROUTINE NRMCON

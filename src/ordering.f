      Subroutine Orderdata(mat, n, k, rowscore)
C23456789012345678901234567890123456789012345678901234567890123456789012
C d = a distance matrix
C n = number of objects (rows of the data files)
C    
C Based on Pierre Legendre Fortran code and Table 9.5 and 9.10
C We compute the principal coordinate of the first axis only.
C Set the precision level for eigenvalue estimation
      Integer mat(n,k)
      Real sumrow(n),sumtot
      real*8 rowscore(n),colscore(n),toler,epsilon
      epsilon=0.000001
      toler=  0.000001
      if(n.gt.1000) then
         epsilon=0.00001
         toler=  0.00001
      endif
      if(n.gt.10000) then
         epsilon=0.0001
         toler=  0.0001
      endif
C Centre the distance matrix (Gower's method)
      call Centre(mat, n, k, sumrow, sumtot)
      call TWWS(mat,n,k,sumrow,sumtot,rowscore,colscore,toler,epsilon)
      end

      Subroutine Centre(mat, n, k, sumrow, sumtot)
      Integer mat(n,k)
      Real d
      Real sumrow(n),sumtot
      do i=1,n
         sumrow(i)=0.0
         enddo
      do i=1,(n-1)
         do j=(i+1),n
            call SM(mat, n, k, i, j, d)
            d = -0.5*d**2
            sumrow(i)=sumrow(i) + d
            sumrow(j)=sumrow(j) + d
            enddo
         enddo
C     
      fln=1.0/float(n)
      sumtot=0.
      do i=1,n
         sumtot=sumtot+sumrow(i)
         sumrow(i)=sumrow(i)*fln
         enddo
      sumtot=sumtot/float(n**2)
C      do 16 i=1,n
C      do 16 j=1,n
C   16 d(i,j)=d(i,j)-sumrow(i)-sumrow(j)+sumtot
      end

      Subroutine SM(mat, n, k, i, j, d)
C Compute a simple matching coefficient from a table of K-means results (integers).
C The 'n' rows are the objects; the 'k' columns are the partitions.
      Integer mat(n,k)
      Real d
C
      a = 0.0
      do kk=1,k
         if(mat(i,kk).eq.mat(j,kk)) a=a+1
         enddo
      d = 1.0 - (a/float(k))
      end

      Subroutine TWWS(mat,n,k,sumrow,sumtot,rowscore,colscore,
     +                toler,epsilon)
      Integer n, niter
      Real sumrow(n),sumtot,d
      Real*8 rowscore(n),colscore(n),epsilon,oldS,newS,toler,
     +       oldrowsc(n)
      niter=1000
C      Step 2: Take the column order as arbitrary initial site scores
      do 4 i=1,n
    4 rowscore(i)=dfloat(i)
      oldS=0.
C      Iterations starting
      do 20 it=1,niter
C      Step 3: Calculate new variable scores (equation 5.8, p. 119)
      do 6 i=1,n
    6 colscore(i)=rowscore(i)  
C      Step 4: Calculate new site scores (equation 5.9, p. 122)
      do 8 i=1,n
      rowscore(i)=0.
      do 8 j=1,n
         call SM(mat, n, k, i, j, d)
         d = -0.5*d**2
         d = d-sumrow(i)-sumrow(j)+sumtot
    8    rowscore(i)=rowscore(i)+d*colscore(j)
C      Step 6: Normalize the site scores
      call NormTWWS(rowscore,n,newS)
      if(newS.lt.epsilon) then
         write(*,103) 0
         goto 52
         endif
C When convergence has been attained, check sign of eigenvalue.
C If ALL rowscores have changed sign during the last iteration, 
C this is a negative eigenvalue.
C      write(*,*) 'oldS-newS', oldS-newS
      if(dabs(oldS-newS).le.toler) then
C        write(*,105) toler, it 
        goto 22
        endif
      do 18 i=1,n
   18 oldrowsc(i)=rowscore(i)
      oldS=newS
   20 continue
C      End of iterations for estimating eigenvalues and eigenvectors
C      write(*,101) k
   22 continue
C      Save Eigenvalues, % variance, cumulative % variance
C      write(*,*) 'Eigenvalue =', newS
C      End of main loop on axes
   52 continue
C
C      Normalize the principal coordinates to variance = eigenvalue
      do 60 i=1,n
   60 rowscore(i)=rowscore(i)*dsqrt(newS)
C      write(*,*) rowscore   
C  101 format(' Convergence not reached for axis:',i3/
C     +       ' Increase NITER or lower TOLER')
C  102 format(' N. iterations to reach convergence for axis',i3,' =',i4)
  103 format(' There are',i4,' eigenvalues different from 0')
C  104 format(' Eigenvector',i3,' is complex [multiply values*Sqrt(-1)]')
C  105 format(" Tolerance is: ", F12.8, "  NIter is: ", i4)

C  150 format(1x,1hR,2i3,6f10.5)
C  151 format(1x,1hC,2i3,10x,6f10.5)
      end


      Subroutine NormTWWS(rowscore,n,newS)
      Integer n
      Real*8 rowscore(n),s2,newS
C      Normalization for two-way weighted summation algorithm for PCA
C      (ter Braak 1987: 123)
C      On output, vector 'rowscore' has length = 1.
C      S = eigenvalue*(p-1) = contraction of vector rowscore in final
C      iteration
C
      s2=0.0
      do 10 i=1,n
   10 s2=s2+rowscore(i)**2
      newS=dsqrt(s2)
      do 20 i=1,n
   20 rowscore(i)=rowscore(i)/newS
      return
      end

c Decorana: the classic goodie by Mark O. Hill
c Ported to R by Jari Oksanen, September 2001:
c * Kept only the numerical engine: subroutines eigy, cutup and
c   yxmult, and the subroutines eigy calls.
c * Changed to lowercase so that it looks prettier.
c * Changed to double precision throughout, with 
c   'implicit double precision (a-h,o-z)' for automatic variables.
c * Changed 'ifix' to 'int', since g77 doesn't know 'ifix' for double.
c * Removed all write statements (they don't show in GUI term --
c   else it would make sense to add a `trace' parameter.)
c
c ORIGINAL INTRODUCTORY COMMENTS:
c
c Cornell Ecology Program Decorana - written by M.O. Hill, July 1979
c Source code and accompanying documentation are available from
c Hugh G. Gauch, Jr., Ecology and Systematics, Cornell University,
c Ithaca, New York 14850.
c
c Performs detrended correspondence analysis;  also will do reciprocal
c averaging as a special case.
c
c **** further modified by Dr Peter R. Minchin May-June 1997
c      - changed tolerance to 0.000005 and iteration limit to 999 in
c          eigy - strict settings of oksanen & minchin 1997
c      - corrected the order-dependent bug in smooth
c      - added parameter statements for maximum dimensions
c      - changed max. dimensions to 5000, 5000, 330000
c      - updated handling of file names retrieved from command line
c      - moved character data to character variables
c      - now accepts relaxed ccf format, with maximum number of pairs
c          per record on line 3
c see: Oksanen, J. & Minchin, P.R. 1997. Instability of ordination
c        results under changes in input data order: explanations and
c        remedies. Journal of Vegetation Science 8: 447-454.
c
c
      subroutine cutup(x,ix,mi,mk)
c takes a vector x and cuts up into (mk-4) segments, putting a
c segmented version of the vector in ix.
      implicit double precision (a-h,o-z)
      double precision x(mi)
      integer ix(mi)
      mmk=mk-4
      maxk=mk-2
      call xmaxmi(x,axmax,axmin,mi)
      axbit=(axmax-axmin)/float(mmk)
      do 10 i=1,mi
C---      iax=ifix((x(i)-axmin)/axbit)+3
      iax=int((x(i)-axmin)/axbit)+3
      if(iax.lt.3) iax=3
      if(iax.gt.maxk) iax=maxk
   10 ix(i)=iax
      return
      end
c
      subroutine trans(y,yy,
     1x,neig,ira,aidot,xeig1,xeig2,xeig3,ix1,ix2,ix3,
     2mi,mk,n,nid,ibegin,iend,idat,qidat)
c this subroutine is the crux of the whole program, in that it
c takes a set of species scores y and iterates to find a new set
c of scores yy.  repeated iteration of this subroutine would lead
c eventually to the correct solution (except that the scores need
c to be divided by the y-totals adotj at each iteration).  the
c calling program eigy is made lengthy by some fancy algebra put
c there to speed up the calculation.  essentially trans is the
c standard reciprocal averaging iteration with either detrending
c with respect to previously derived axes (in the case of detrended
c correspondence analysis) or orthogonalization with respect to
c them (in the case of reciprocal averaging).
      implicit double precision (a-h,o-z)
      double precision x(mi),xeig1(mi),xeig2(mi),xeig3(mi)
      double precision y(n),yy(n),aidot(mi),qidat(nid)
      integer ix1(mi),ix2(mi),ix3(mi),idat(nid),ibegin(mi),iend(mi)
      call yxmult(y,x,mi,n,nid,ibegin,iend,idat,qidat)
      do 10 i=1,mi
   10 x(i)=x(i)/aidot(i)
      if(neig.eq.0) goto 200
      if(ira.eq.1) goto 100
      call detrnd(x,aidot,ix1,mi,mk)
      if(neig.eq.1) goto 200
      call detrnd(x,aidot,ix2,mi,mk)
      if(neig.eq.2) goto 90
      call detrnd(x,aidot,ix3,mi,mk)
      call detrnd(x,aidot,ix2,mi,mk)
   90 call detrnd(x,aidot,ix1,mi,mk)
      goto 200
  100 a1=0.0
      do 110 i=1,mi
  110 a1=a1+aidot(i)*x(i)*xeig1(i)
      do 120 i=1,mi
  120 x(i)=x(i)-a1*xeig1(i)
      if(neig.eq.1) goto 200
      a2=0.0
      do 130 i=1,mi
  130 a2=a2+aidot(i)*x(i)*xeig2(i)
      do 140 i=1,mi
  140 x(i)=x(i)-a2*xeig2(i)
      if(neig.eq.2) goto 200
      a3=0.0
      do 150 i=1,mi
  150 a3=a3+aidot(i)*x(i)*xeig3(i)
      do 160 i=1,mi
  160 x(i)=x(i)-a3*xeig3(i)
  200 call xymult(x,yy,mi,n,nid,ibegin,iend,idat,qidat)
      return
      end
c
      subroutine detrnd(x,aidot,ix,mi,mk)
      implicit double precision (a-h,o-z)
      double precision x(mi),z(50),zn(50),zbar(50),aidot(mi)
      integer ix(mi)
c starts with a vector x and detrends with respect to groups defined
c by ix.  detrending is in blocks of 3 units at a time, and the
c result calculated is the average of the 3 possible results that
c can be obtained, corresponding to 3 possible starting positions
c for the blocks of 3.
      do 10 k=1,mk
      z(k)=0.0
   10 zn(k)=0.0
      do 20 i=1,mi
      k=ix(i)
      z(k)=z(k)+x(i)*aidot(i)
   20 zn(k)=zn(k)+aidot(i)
      mmk=mk-1
      do 30 k=2,mmk
   30 zbar(k)=(z(k-1)+z(k)+z(k+1))/(zn(k-1)+zn(k)+zn(k+1)+1.0e-12)
      mmk=mmk-1
      do 35 k=3,mmk
   35 z(k)=(zbar(k-1)+zbar(k)+zbar(k+1))/3.0
      do 40 i=1,mi
      k=ix(i)
   40 x(i)=x(i)-z(k)
      return
      end
c
      subroutine eigy(x,y,eig,neig,ira,iresc,short,
     1mi,mk,n,nid,ibegin,iend,idat,qidat,y2,y3,y4,y5,
     2xeig1,xeig2,xeig3,ix1,ix2,ix3,aidot,adotj)
c extracts an eigenvector y with eigenvalue eig.  the algebra
c is a little complicated, but consists essentially of repre-
c senting the transformation (subroutine trans) approximately
c by a tridiagonal 4x4 matrix.  the eigenproblem for the
c tridiagonal matrix is solved and this solution is plugged
c back in to obtain a new trial vector.
c after getting the eigenvector, the scores may be rescaled
c (subroutine strtch).
      implicit double precision (a-h,o-z)
      double precision x(mi),y(n),y2(n),y3(n),y4(n),y5(n)
      double precision xeig1(mi),xeig2(mi),xeig3(mi),aidot(mi),adotj(n)
      double precision qidat(nid)
      integer ibegin(mi),iend(mi),idat(nid),ix1(mi),ix2(mi),ix3(mi)
c string to print R warnings: this must be long enough to fit format
c statement 1012
      character*64 warning
      tot=0.0
      do 10 j=1,n
      tot=tot+adotj(j)
   10 y(j)=float(j)
      y(1)=1.1
c---tolerance reduced by p.minchin jan 1997
c      tol=0.0001
      tol=0.000005
      call trans(y,y,
     1x,neig,ira,aidot,xeig1,xeig2,xeig3,ix1,ix2,ix3,
     2mi,mk,n,nid,ibegin,iend,idat,qidat)
      icount=0
   20 a=0.0
      do 30 j=1,n
   30 a=a+y(j)*adotj(j)
      a=a/tot
      ex=0.0
      do 40 j=1,n
      ay=y(j)-a
      ex=ex+ay*ay*adotj(j)
   40 y(j)=ay
      ex=sqrt(ex)
      do 50 j=1,n
   50 y(j)=y(j)/ex
      call trans(y,y2,
     1x,neig,ira,aidot,xeig1,xeig2,xeig3,ix1,ix2,ix3,
     2mi,mk,n,nid,ibegin,iend,idat,qidat)
      a=0.0
      a11=0.0
      a12=0.0
      a22=0.0
      a23=0.0
      a33=0.0
      a34=0.0
      a44=0.0
      do 60 j=1,n
      ay=y2(j)
      y2(j)=ay/adotj(j)
      a=a+ay
   60 a11=a11+ay*y(j)
      a=a/tot
      do 70 j=1,n
      ay=y2(j)-(a+a11*y(j))
      a12=a12+ay*ay*adotj(j)
   70 y2(j)=ay
      a12=sqrt(a12)
      do 80 j=1,n
   80 y2(j)=y2(j)/a12
c---------removed write statements--------------------------
c      if(icount.eq.0) write(*,1000)
c 1000 format(/1x)
c      write(*,1011) a12,icount
c 1011 format(1x,'residual',f10.6,'       at iteration',i3)
      if(a12.lt.tol) goto 200
c---maximum iteration limit increased by P.Minchin jan 1997
c      if(icount.gt.9) goto 200
      if(icount.gt.999) goto 200
      icount=icount+1
      call trans(y2,y3,
     1x,neig,ira,aidot,xeig1,xeig2,xeig3,ix1,ix2,ix3,
     2mi,mk,n,nid,ibegin,iend,idat,qidat)
      a=0.0
      b13=0.0
      do 90 j=1,n
      ay=y3(j)
      y3(j)=ay/adotj(j)
      a=a+ay
      a22=a22+ay*y2(j)
   90 b13=b13+ay*y(j)
      a=a/tot
      do 100 j=1,n
      ay=y3(j)-(a+a22*y2(j)+b13*y(j))
      a23=a23+ay*ay*adotj(j)
  100 y3(j)=ay
      a23=sqrt(a23)
      if(a23.gt.tol) goto 105
      a23=0.0
      goto 160
  105 continue
      do 110 j=1,n
  110 y3(j)=y3(j)/a23
      call trans(y3,y4,
     1x,neig,ira,aidot,xeig1,xeig2,xeig3,ix1,ix2,ix3,
     2mi,mk,n,nid,ibegin,iend,idat,qidat)
      a=0.0
      b14=0.0
      b24=0.0
      do 120 j=1,n
      ay=y4(j)
      y4(j)=y4(j)/adotj(j)
      a=a+ay
      a33=a33+ay*y3(j)
      b14=b14+ay*y(j)
  120 b24=b24+ay*y2(j)
      a=a/tot
      do 130 j=1,n
      ay=y4(j)-(a+a33*y3(j)+b14*y(j)+b24*y2(j))
      a34=a34+ay*ay*adotj(j)
  130 y4(j)=ay
      a34=sqrt(a34)
      if(a34.gt.tol) goto 135
      a34=0.0
      goto 160
  135 continue
      do 140 j=1,n
  140 y4(j)=y4(j)/a34
      call trans(y4,y5,
     1x,neig,ira,aidot,xeig1,xeig2,xeig3,ix1,ix2,ix3,
     2mi,mk,n,nid,ibegin,iend,idat,qidat)
      do 150 j=1,n
  150 a44=a44+y4(j)*y5(j)
c we now have the tridiagonal representation of trans.  solve
c eigenproblem for tridiagonal matrix.
  160 ax1=1.0
      ax2=0.1
      ax3=0.01
      ax4=0.001
      do 170 itimes=1,100
      axx1=a11*ax1+a12*ax2
      axx2=a12*ax1+a22*ax2+a23*ax3
      axx3=a23*ax2+a33*ax3+a34*ax4
      axx4=a34*ax3+a44*ax4
      ax1=a11*axx1+a12*axx2
      ax2=a12*axx1+a22*axx2+a23*axx3
      ax3=a23*axx2+a33*axx3+a34*axx4
      ax4=a34*axx3+a44*axx4
      ex=sqrt(ax1**2+ax2**2+ax3**2+ax4**2)
      ax1=ax1/ex
      ax2=ax2/ex
      ax3=ax3/ex
      ax4=ax4/ex
      if(itimes.ne.(itimes/5)*5) goto 170
      exx=sqrt(ex)
      resi=sqrt((ax1-axx1/exx)**2+(ax2-axx2/exx)**2+
     1(ax3-axx3/exx)**2+(ax4-axx4/exx)**2)
      if(resi.lt.tol*0.05) goto 180
  170 continue
  180 continue
      do 190 j=1,n
  190 y(j)=ax1*y(j)+ax2*y2(j)+ax3*y3(j)+ax4*y4(j)
      goto 20
c-----------Removed write statements, added 200 continue--------
 200  continue
c  200 write(*,1010) a11
c 1010 format(1x,'eigenvalue',f10.5)
c      if(a12.gt.tol) write(*,1012) tol
c 1012 format(1x,'*** beware ***     residual bigger than tolerance',
c     1', which is',f10.6)
c R version of the above warning. You must change the length of
c character*n warning definition if you change the warning text
      if (a12 .gt. tol) then
         write(warning, 1012) a12, tol, neig+1
 1012    format("residual", f10.7, " bigger than tolerance", f10.7, 
     1 " for axis ", i1)
         call rwarn(warning)
      end if
c we calculate x from y, and set x to unit length if reciprocal
c averaging option is in force (ira=1)
      call xmaxmi(y,aymax,aymin,n)
      sign=1.0
      if(-aymin.gt.aymax) sign=-1.0
      do 210 j=1,n
  210 y(j)=y(j)*sign
      call yxmult(y,x,mi,n,nid,ibegin,iend,idat,qidat)
      do 220 i=1,mi
  220 x(i)=x(i)/aidot(i)
      if(iresc.eq.0) goto 225
      if(a11.gt.0.999) goto 225
      do 223 i=1,iresc
      monit=0
      if(i.eq.1.or.i.eq.iresc) monit=1
      call strtch(x,y,short,monit,
     1mi,n,nid,aidot,ibegin,iend,idat,qidat)
  223 continue
      eig=a11
      return
  225 axlong=0.0
      do 230 i=1,mi
  230 axlong=axlong+aidot(i)*x(i)**2
      axlong=sqrt(axlong)
      do 240 i=1,mi
  240 x(i)=x(i)/axlong
      do 250 j=1,n
  250 y(j)=y(j)/axlong
c it remains to scale y to unit within-sample standard deviation
      sumsq=0.0
      do 260 i=1,mi
      id1=ibegin(i)
      id2=iend(i)
      ax=x(i)
      do 255 id=id1,id2
      j=idat(id)
  255 sumsq=sumsq+qidat(id)*(ax-y(j))**2
  260 continue
      sd=sqrt(sumsq/tot)
      if(a11.le.0.999) goto 265
      sd=aymax/axlong
      sd1=-aymin/axlong
      if(sd1.gt.sd) sd=sd1
  265 continue
      do 270 j=1,n
  270 y(j)=y(j)/sd
      eig=a11
      return
      end
c
      subroutine xymult(x,y,mi,n,nid,ibegin,iend,idat,qidat)
c starts with vector x and forms matrix product y=ax
      implicit double precision (a-h,o-z)
      double precision x(mi),y(n),qidat(nid)
      integer ibegin(mi),iend(mi),idat(nid)
      do 10 j=1,n
   10 y(j)=0.0
      do 30 i=1,mi
      id1=ibegin(i)
      id2=iend(i)
      ax=x(i)
      do 20 id=id1,id2
      j=idat(id)
   20 y(j)=y(j)+ax*qidat(id)
   30 continue
      return
      end
c
      subroutine yxmult(y,x,mi,n,nid,ibegin,iend,idat,qidat)
c starts with vector y and forms matrix product x=ay
      implicit double precision (a-h,o-z)
      double precision x(mi),y(n),qidat(nid)
      integer ibegin(mi),iend(mi),idat(nid)
      do 20 i=1,mi
      id1=ibegin(i)
      id2=iend(i)
      ax=0.0
      do 10 id=id1,id2
      j=idat(id)
   10 ax=ax+y(j)*qidat(id)
   20 x(i)=ax
      return
      end
c
      subroutine xmaxmi(x,axmax,axmin,m)
c forms maximum and minimum of x(m)
      implicit double precision (a-h,o-z)
      double precision x(m)
      axmax=-1.0e10
      axmin=-axmax
      do 10 i=1,m
      ax=x(i)
      if(ax.gt.axmax) axmax=ax
      if(ax.lt.axmin) axmin=ax
   10 continue
      return
      end
c
c
      subroutine strtch(x,y,short,monit,
     1mi,n,nid,aidot,ibegin,iend,idat,qidat)
c takes an axis (x,y) and scales to unit mean square dev of species
c scores per sample.  an attempt is made for longer axes (l > short)
c to squeeze them in and out so that they have the right mean square
c deviation all the way along the axis and not only on average.
      implicit double precision (a-h,o-z)
      double precision x(mi),y(n),aidot(mi),qidat(nid)
      double precision zn(50),zv(50)
      integer ibegin(mi),iend(mi),idat(nid)
c Use monit so that gfortran -Wall does not warn for it being unused
      monit = 0
c---common block added by p.minchin feb 1988
c      common /lunits/ iuinp1,iuinp2,iuout1,iuout2,iuout3
      do 200 icount=1,2
      mk=20
      call segmnt(x,y,zn,zv,mi,mk,n,nid,
     1aidot,ibegin,iend,idat,qidat)
      call smooth(zv,mk)
      call smooth(zn,mk)
      zvsum=0.0
      do 50 k=1,mk
   50 zvsum=zvsum+zv(k)/zn(k)
      sd=sqrt(zvsum/float(mk))
c we want mean within-sample square deviation to be 1.0, so we divide
c everything by sd
      along=0.0
      do 60 i=1,mi
      ax=x(i)/sd
      x(i)=ax
   60 if(along.lt.ax) along=ax
c--------Removed write statements----------------------------
c      if(icount.eq.1.and.monit.eq.1) write(*,1000)
c 1000 format(/1x)
c      if(monit.eq.1) write(*,1001) along
c 1001 format(1x,'length of gradient',f10.3)
      do 70 j=1,n
   70 y(j)=y(j)/sd
      if(along.lt.short) return
      if(icount.eq.2) return
c      mk=ifix(along*5.0)+1
      mk=int(along*5.0)+1
      if(mk.lt.10) mk=10
      if(mk.gt.45) mk=45
      call segmnt(x,y,zn,zv,mi,mk,n,nid,aidot,ibegin,iend,idat,qidat)
      call smooth(zv,mk)
      call smooth(zn,mk)
      zvsum=0.0
      do 100 k=1,mk
      azv=1.0/sqrt(0.2/along+zv(k)/zn(k))
      zvsum=zvsum+azv
  100 zv(k)=azv
      do 110 k=1,mk
  110 zv(k)=zv(k)*along/zvsum
c----------Removed write statements-------------------------
c      if(monit.eq.1) write(*,1002) (zv(k),k=1,mk)
c 1002 format(1x,'length of segments',10f6.2)
      az=0.0
      zn(1)=0.0
      do 120 k=1,mk
      az=az+zv(k)
  120 zn(k+1)=az
      axbit=along/float(mk)
      do 130 j=1,n
C      iay=ifix(y(j)/axbit)+1
      iay=int(y(j)/axbit)+1
      if(iay.lt.1) iay=1
      if(iay.gt.mk) iay=mk
  130 y(j)=zn(iay)+zv(iay)*(y(j)/axbit-float(iay-1))
      call yxmult(y,x,mi,n,nid,ibegin,iend,idat,qidat)
      do 140 i=1,mi
  140 x(i)=x(i)/aidot(i)
  200 continue
      return
      end
c
      subroutine smooth(z,mk)
      implicit double precision (a-h,o-z)
      double precision z(mk)
c takes a vector z and does (1,2,1)-smoothing until no blanks left
c and then 2 more iterations of (1,2,1)-smoothing.  if no blanks to
c begin with, then does 3 smoothings, i.e. effectively (1,6,15,20,
c 15,6,1)-smoothing.
      istop=1
      do 20 icount=1,50
      az2=z(1)
      az3=z(2)
      if(az3.eq.0.0) istop=0
      z(1)=0.75*az2+0.25*az3
      do 10 k3=3,mk
      az1=az2
      az2=az3
      az3=z(k3)
c---bug in next line fixed by p.minchin jan 1997
c      if(az3.lt.0.0) istop=0
      if(az3.le.0.0) istop=0
   10 z(k3-1)=0.5*(az2+0.5*(az1+az3))
      z(mk)=0.25*az2+0.75*az3
      istop=istop+1
      if(istop.eq.4) goto 30
   20 continue
   30 return
      end
c
      subroutine segmnt(x,y,zn,zv,mi,mk,n,nid,
     1aidot,ibegin,iend,idat,qidat)
c given an ordination (x,y), calculates numbers and summed mean-square
c deviations in mk segments.  zn(k) is the number of samples in segment
c k;  zv(k) is the summed mean-square deviation.  (we aim to make zv,
c zn as nearly equal as possible.)
      implicit double precision (a-h,o-z)
      double precision x(mi),y(n),zn(mk),zv(mk),aidot(mi),qidat(nid)
      integer ibegin(mi),iend(mi),idat(nid)
      do 10 k=1,mk
      zn(k)=-1.0e-20
   10 zv(k)=-1.0e-20
      call xmaxmi(x,axmax,axmin,mi)
      axbit=(axmax-axmin)/float(mk)
      do 20 i=1,mi
   20 x(i)=x(i)-axmin
      do 30 j=1,n
   30 y(j)=y(j)-axmin
      do 50 i=1,mi
      sqcorr=0.0
      sumsq=2.0e-20
      id1=ibegin(i)
      id2=iend(i)
      ax=x(i)
      do 40 id=id1,id2
      j=idat(id)
      aij=qidat(id)
      sqcorr=sqcorr+aij**2
   40 sumsq=sumsq+aij*(ax-y(j))**2
      sqcorr=sqcorr/aidot(i)**2
      if(sqcorr.gt.0.9999) sqcorr=0.9999
      sumsq=sumsq/aidot(i)
c      k=ifix(ax/axbit)+1
      k=int(ax/axbit)+1
      if(k.gt.mk) k=mk
      if(k.lt.1) k=1
      zv(k)=zv(k)+sumsq
      zn(k)=zn(k)+1.0-sqcorr
   50 continue
      return
      end


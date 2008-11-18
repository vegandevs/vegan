C *************************************************************
C Based on my standard cepin module. This had to be vandalized,
C since R cannot interface with character arrays: It can handle only
C a single character*255 variable. So the routine had to be broken
C into pieces, and assume that the file opened  is still there when 
C called repeatedly. 
C
C I feel really bad for this.
C
C Since the subroutine had to be split, I split it completely:
C Each file type is read in a separate subroutine, so that I can
C get rid off static work arrays.
C
C Jari Oksanen, Sep 16, 2001
C
C *************************************************************

C
C Read the file header.
C The file is left open (that's magic...), and everything
C needed is the format, the `kind' of the file, the number
C of items in a record (`nitem'), and number of records for 
C `kind=1'. The `kind' is interpreted from number of I's in `fmt'. 
C

      subroutine cephead(cepfile, kind, nitem, nst, fmt)

      character*255 cepfile, fmt
      integer kind, nitem, nst

      character*80 title
      integer ni, nbrac

      ni=0
      nbrac=0

      open (unit=1,file=cepfile,status='old')
      read (1,101) title
 101  format (a80)
      read (1,101) fmt
      do i=1,80
         if (fmt(i:i) .eq. 'I' .or. fmt(i:i) .eq. 'i') ni=ni+1
      end do
      kind = ni+1
      if (kind .eq. 1) then
         read(1, *) nitem, nst
         return
      endif
      nbrac=0
      do i=1,80
         if (fmt(i:i) .eq. '(' ) nbrac=nbrac+1
         if (fmt(i:i) .eq. ')' ) then
            nbrac=nbrac-1
            if (nbrac .le. 0) goto 6
         end if
         goto 8
 6       ii=i+2
         if (ii .le. 69) read(fmt,'(t69,i2)') nitem
         if (nitem .eq. 0) then
            read (1,*) nitem
         else
            fmt(ii:80)=' '
         endif
 8       continue
      enddo
      return
      end
      
c
c Open CEP format (kind=2).
c Zeros are skipped, negative entries stored
c Stop when negative site (row) index found.
c
      subroutine cepopen(fmt, nitem, maxdat, nsp, nst, idplot, idspec, 
     X abund, work, ier)

      character*255 fmt
      integer nitem, nsp, nst
      integer idplot(maxdat), idspec(maxdat)
      double precision abund(maxdat)
 
      integer id, ii, j, ier
      double precision work(nitem)

      nsp = nitem
      nst = 0
      id = 0
      ier = 99
      
 20   read (1, fmt) ii, (work(j),j=1,nitem)  
      if (ii .le. 0) then
         ier = 0
         return
      endif
      if (ii .gt. nst) then
         nst=ii
      end if
      do j=1,nitem
         if (work(j) .ne. 0) then
            id = id+1
            if (id .gt. maxdat) then
               ier = 1
               return
            endif
            idplot(id) = ii
            idspec(id) = j
            abund(id) = work(j)
         endif
      end do
      goto 20
      
      end

C
C Condensed CEP format (kind=3)
C ALL entries are stored in condensed format (except zeros)
C

      subroutine cepcond(fmt, nitem, maxdat, nsp, nst, idplot, idspec, 
     X abund, work, item, ier)

      character*255 fmt
      integer nitem, nsp, nst
      integer idplot(maxdat), idspec(maxdat)
      double precision abund(maxdat)
 
      integer id, ii, j, ier
      integer item(nitem)
      double precision work(nitem)

      nsp = 0
      nst = 0
      id = 0
      ier = 99

 40   read (1,fmt) ii,(item(j),work(j),j=1,nitem) 
      if (ii .le. 0) then
         ier = 0
         return
      endif
      if (ii .gt. nst) then
         nst = ii
      endif
      do j = 1,nitem
         if (item(j) .gt. 0 .and. work(j) .ne. 0.0) then
            id = id+1
            if (id .gt. maxdat) then
               ier = 1
               return
            endif
            idplot(id) = ii
            if (item(j) .gt. nsp) then
               nsp = item(j)
            endif
            idspec(id) = item(j)
            abund(id) = work(j)
         endif
      enddo
      goto 40

      end

c      
c Free CEP format (kind=1) -- or error!
c Get everything but zeros.
c

      subroutine cepfree(nitem, maxdat, nsp, nst, idplot, idspec, 
     X abund, work, ier)

      integer nitem, nsp, nst, maxdat
      integer idplot(maxdat), idspec(maxdat)
      double precision abund(maxdat)
 
      integer id, j, i, ier
      double precision work(nitem)

      nsp = nitem
      nst = nst
      id = 0
      ier = 99

      do i=1,nst
         read(1,*) (work(j),j=1,nsp)
         do j=1,nsp
            if(work(j)  .ne. 0) then
               id=id+1
               if (id .gt. maxdat) then
                  ier = 1
                  return
               endif
               abund(id) = work(j)
               idspec(id) = j
               idplot(id) = i
            endif
         enddo
      enddo
      
      ier = 0

      return
      end

c
c     R can handle only a single character*255: 
c     Names cannot be read in one time, but they must be
c     handled line by line.
c

      subroutine cepnames(entname)
      character*255 entname
      read (1,1100, end=666) entname  
 1100 format (a80)
 666  continue
      return
      end

c
c Stupid, but needed
c

      subroutine cepclose()
      close(1)
      return
      end


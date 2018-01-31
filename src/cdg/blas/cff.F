c======================================================================================
c Name         : cff
c Date         : 1/1/18 (DLR)
c Copyright    : Copyright 2018-2018. Colorado State University. All rights reserved
c Description  : Cache-friendly optimized Fortran linear algebra routines
c======================================================================================




c************************************************************************************
c************************************************************************************
c METHOD     : dmxm
c DESCRIPTION: mxm driver routine that that calls methods explicitly unrolled
c ARGUMENTS  : 
c RETURNS    : 
c************************************************************************************
      subroutine dmxm(C, A, nai, naj, B, nbi, nbj, isz)
      real*8  C(nai,nbj), A(nai,naj), B(nbi,nbj)
      real*8  alpha, beta

      if      ( nbi .eq.  1 ) then
        call dmxm1 (A, nai, B, nbj, C, nbj)
      else if ( nbi .eq.  2 ) then
        call dmxm2 (A, nai, B, nbj, C, nbj)
      else if ( nbi .eq.  3 ) then
        call dmxm3 (A, nai, B, nbj, C, nbj)
      else if ( nbi .eq.  4 ) then
        call dmxm4 (A, nai, B, nbj, C, nbj)
      else if ( nbi .eq.  5 ) then
        call dmxm5 (A, nai, B, nbj, C, nbj)
      else if ( nbi .eq.  6 ) then
        call dmxm6 (A, nai, B, nbj, C, nbj)
      else if ( nbi .eq.  7 ) then
        call dmxm7 (A, nai, B, nbj, C, nbj)
      else if ( nbi .eq.  8 ) then
        call dmxm8 (A, nai, B, nbj, C, nbj)
      else if ( nbi .eq.  9 ) then
        call dmxm9 (A, nai, B, nbj, C, nbj)
      else if ( nbi .eq. 10 ) then
        call dmxm10(A, nai, B, nbj, C, nbj)
      else if ( nbi .eq. 11 ) then
        call dmxm11(A, nai, B, nbj, C, nbj)
      else if ( nbi .eq. 12 ) then
        call dmxm12(A, nai, B, nbj, C, nbj)
      else if ( nbi .eq. 13 ) then
        call dmxm13(A, nai, B, nbj, C, nbj)
      else if ( nbi .eq. 14 ) then
        call dmxm14(A, nai, B, nbj, C, nbj)
      else if ( nbi .eq. 15 ) then
        call dmxm15(A, nai, B, nbj, C, nbj)
      else if ( nbi .eq. 16 ) then
        call dmxm16(A, nai, B, nbj, C, nbj)
      else if ( nbi .eq. 17 ) then
        call dmxm17(A, nai, B, nbj, C, nbj)
      else if ( nbi .eq. 18 ) then
        call dmxm18(A, nai, B, nbj, C, nbj)
      else if ( nbi .eq. 19 ) then
        call dmxm19(A, nai, B, nbj, C, nbj)
      else if ( nbi .eq. 20 ) then
        call dmxm20(A, nai, B, nbj, C, nbj)
      else if ( nbi .eq. 21 ) then
        call dmxm21(A, nai, B, nbj, C, nbj)
      else if ( nbi .eq. 22 ) then
        call dmxm22(A, nai, B, nbj, C, nbj)
      else if ( nbi .eq. 23 ) then
        call dmxm23(A, nai, B, nbj, C, nbj)
      else if ( nbi .eq. 24 ) then
        call dmxm24(A, nai, B, nbj, C, nbj)
      else
#if defined(_GBLAS)
        call dmxmcf(C, A, nai, naj, B, nbi, nbj, isz )
#else
        alpha = 1.0
        beta  = 0.0
        call dgemm("N","N", nai, nbj, alpha, A, nai, 
     1             B, nbi, beta, C, nai)
#endif
      endif


      return
c
c  end of subroutine dmxm
c
      end


c************************************************************************************
c************************************************************************************
c METHOD     : dmxmcf
c DESCRIPTION: Cache-friendly mxm routine that does not unroll by hand
c ARGUMENTS  : 
c RETURNS    : 
c************************************************************************************
      subroutine dmxmcf(C, A, nai, naj, B, nbi, nbj, isz)
      real*8  C(nai,nbj), A(nai,naj), B(nbi,nbj)
      real*8  zero
      integer nai, naj, nbi, nbj, isz
      integer i, ii, j, jj, k, kk

c  Initialize C to 0:
      zero = 0.0
      call matset(C, nai, nbj, zero, isz)

c  Do matrix-matrix multiply using cache-blocking:
      do ii = 1, nai, isz
        do jj = 1, nbj, isz
          do kk = 1, nbi, isz
  
            do i = ii, min(nai,ii+isz-1)
              do j = jj, min(nbj,jj+isz-1)
                do k = kk, min(nbi,kk+isz-1)
                  C(i,j) = C(i,j) + A(i,k)*B(k,j)
                enddo
              enddo
            enddo
        
          enddo
        enddo
      enddo


      return
c
c  end of subroutine dmxmfc
c
      end

c************************************************************************************
c************************************************************************************
c METHOD     : dmxv
c DESCRIPTION:
c ARGUMENTS  : 
c RETURNS    : 
c************************************************************************************
      subroutine dmxv(y, ny, A, x, nx, isz)
      real*8  y(*), A(ny,*), x(*)
      integer ny, nx, isz
      integer i, ii, j, jj

      do ii = 1, ny, isz
        do i = ii, min(ny,ii+isz-1)
           y(i) = 0.0
        enddo
      enddo


c  Do matrix-vector multiply:
      do ii = 1, ny, isz
        do jj = 1, nx, isz

          do i = ii, min(ny,ii+isz-1)
            do j = jj, min(nx,jj+isz-1)
              y(i) = y(i) + A(i,j)*x(j)
            enddo
          enddo

        enddo
      enddo

      return
c
c  end of subroutine dmxv
c
      end

c************************************************************************************
c************************************************************************************
c METHOD     : dmxDm  
c DESCRIPTION: computes C = A * Diag(B). Dimensions of C assumed to be 
c              those of A
c ARGUMENTS  : 
c RETURNS    : 
c************************************************************************************
      subroutine dmxDm(C, A, nai, naj, b, nb, isz)
      real*8  C(nai,naj), A(nai,naj), b(*)
      integer nai, naj, nb, isz

      integer i, ii, j, jj

c  Do matrix-Diag(matrix) multiply:
      do ii = 1, nai, isz 
        do jj = 1, naj, isz 

          do i = ii, min(nai,ii+isz-1)
            do j = jj, min(naj,jj+isz-1)
              C(i,j) = A(i,j)*b(j)
            enddo
          enddo

        enddo
      enddo   

      return  
c
c  end of subroutine dmxDm
c
      end     


c************************************************************************************
c************************************************************************************
c METHOD     : dDmxm  
c DESCRIPTION: computes C = Diag(A)*B. Dimensions of C assumed to be 
c              those of B
c ARGUMENTS  : 
c RETURNS    : 
c************************************************************************************
      subroutine dDmxm(C, a, na, B, nbi, nbj, isz)
      real*8  C(nbi,nbj), B(nbi,nbj), a(na)
      integer nbi, nbj, na, isz

      integer i, ii, j, jj

c  Do Diag(matrix)-matrix multiply:
      do ii = 1, nbi, isz 
        do jj = 1, nbj, isz 

          do i = ii, min(nbi,ii+isz-1)
            do j = jj, min(nbj,jj+isz-1)
              C(i,j) = a(i)*B(i,j)
            enddo
          enddo

        enddo
      enddo   

      return  
c
c  end of subroutine dDmxm
c
      end     


c************************************************************************************
c************************************************************************************
c METHOD     : daApbB  
c DESCRIPTION: computes C = a A + b B. Dimensions of A, B, C assumed to be equal
c              (nxm), and a, and b are double constants
c ARGUMENTS  : 
c RETURNS    : 
c************************************************************************************
      subroutine daApbB(C, A, B, n, m, ca, cb, isz)
      real*8  A(n,m), B(n,m), C(n,m), ca, cb
      integer n, m, isz

      integer i, ii, j, jj

c  Do sum:
      do ii = 1, m, isz 
        do jj = 1, n, isz 

          do i = ii, min(m,ii+isz-1)
            do j = jj, min(n,jj+isz-1)
              C(i,j) = ca*A(i,j) + cb*B(i,j)
            enddo
          enddo

        enddo
      enddo   

      return  
c
c  end of subroutine daApbB
c
      end     


c************************************************************************************
c************************************************************************************
c METHOD     : dzaxpby
c DESCRIPTION:
c ARGUMENTS  : 
c RETURNS    : 
c************************************************************************************
      subroutine dzaxpby(z, x, a, y, b, nxy, isz)
      real*8  z(nxy), x(nxy), y(nxy), a, b
      real*8  prod
      integer nxy, isz
      integer i, ii

      isz = isz

c  Do vector sums:
#if 0
      do ii = 1, nxy, isz
        do i = ii, min(nxy,ii+isz-1)
          prod = b*y(i) 
          z(i) = a*x(i) + prod
        enddo
      enddo
#endif

      do i = 1, nxy
        z(i) = a*x(i) + b*y(i)
      enddo

      return
c
c  end of subroutine dzaxpby
c
      end


c************************************************************************************
c************************************************************************************
c METHOD     : dxaxpby
c DESCRIPTION: x(i) = a*x(i)+b*y(i) 
c ARGUMENTS  : 
c RETURNS    : 
c************************************************************************************
      subroutine dxaxpby(x, a, y, b, nxy, isz)
      real*8  x(nxy), y(nxy), a, b
      real*8  prod
      integer nxy, isz
      integer i, ii

c  Do vector sums:
#if 0
      do ii = 1, nxy, isz
        do i = ii, min(nxy,ii+isz-1)
          prod = b*y(i) + 0.0
          x(i) = a*x(i) + prod
        enddo
      enddo
#endif
      do i = 1, nxy
        x(i) = a*x(i) + b*y(i)
      enddo

      return
c
c  end of subroutine dxaxpby
c
      end


c************************************************************************************
c************************************************************************************
c METHOD     : dzvxvpt
c DESCRIPTION: z(i) = x(i)*y(i) 
c ARGUMENTS  : 
c RETURNS    : 
c************************************************************************************
      subroutine dzvxvpt(z, x, y, nxy, isz)
      real*8  z(nxy), x(nxy), y(nxy)
      integer nxy, isz
      integer i, ii

c  Do point product:
#if 0
      do ii = 1, nxy, isz
        do i = ii, min(nxy,ii+isz-1)
          z(i) = x(i) * y(i) 
        enddo
      enddo
#endif

      do i = 1, nxy
        z(i) = x(i) * y(i) 
      enddo

      return
c
c  end of subroutine dzvxvpt
c
      end


c************************************************************************************
c************************************************************************************
c METHOD     : dvvxvpt
c DESCRIPTION: x(i) = x(i)*y(i) 
c ARGUMENTS  : 
c RETURNS    : 
c************************************************************************************
      subroutine dvvxvpt(x, y, nxy, isz)
      real*8  x(nxy), y(nxy)
      integer nxy, isz
      integer i, ii

c  Do point product:
#if 0
      do ii = 1, nxy, isz

        do i = ii, min(nxy,ii+isz-1)
          x(i) = x(i) * y(i) 
        enddo

      enddo
#endif

      do i = 1, nxy
        x(i) = x(i) * y(i) 
      enddo

      return
c
c  end of subroutine dvvxvpt
c
      end


c************************************************************************************
c************************************************************************************
c METHOD     : dvvxvptpv
c DESCRIPTION: x(i) = x(i)*y(i) + cz*z(i)
c ARGUMENTS  : 
c RETURNS    : 
c************************************************************************************
      subroutine dvvxvptpv(x, y, z, cz, nxy, isz)
      real*8  x(nxy), y(nxy), z(nxy), cz
      real*8  prod
      integer nxy, isz
      integer i, ii

c  Do point product & sum:
      do ii = 1, nxy, isz
        do i = ii, min(nxy,ii+isz-1)
          prod = cz*z(i) 
          x(i) = x(i) * y(i) + prod
        enddo
      enddo

      return
c
c  end of subroutine dvvxvptpv
c
      end


c************************************************************************************
c************************************************************************************
c METHOD     : dvpvxvpt
c DESCRIPTION: z(i) = x(i)*y(i) + cz*z(i)
c ARGUMENTS  : 
c RETURNS    : 
c************************************************************************************
      subroutine dvpvxvpt(x, y, z, cz, nxy, isz)
      real*8  x(nxy), y(nxy), z(nxy), cz
      real*8  prod
      integer nxy, isz
      integer i, ii

c  Do point product & sum:
#if 0
      do ii = 1, nxy, isz
        do i = ii, min(nxy,ii+isz-1)
          prod = cz*z(i) 
          z(i) = x(i) * y(i) + prod
        enddo
      enddo
#endif
      do i = 1, nxy
        z(i) = x(i) * y(i) + cz*z(i)
      enddo

      return
c
c  end of subroutine dvpvxvpt
c
      end


c************************************************************************************
c************************************************************************************
c METHOD     : ddotg
c DESCRIPTION:
c ARGUMENTS  : 
c RETURNS    : 
c************************************************************************************
      subroutine ddotg(rdot, x, y, nxy, isz)
      real*8  rdot, x(nxy), y(nxy)
      real*8  e, s, yy, z
      integer nxy, isz
      integer i, ii

c  Do vector dot product:
      rdot = 0.0
      e    = 0.0
      do ii = 1, nxy, isz
        do i = ii, min(nxy,ii+isz-1)
          yy   = x(i)*y(i) - e
          z    = rdot + yy
          e    = (z - rdot) - yy
          rdot = z
        enddo
      enddo

      return
c
c  end of subroutine ddotg
c
      end


c************************************************************************************
c************************************************************************************
c METHOD     : dcopy
c DESCRIPTION: copies y into x
c ARGUMENTS  : 
c RETURNS    : 
c************************************************************************************
      subroutine dcopy(x, y, nxy, isz)
      real*8  x(nxy), y(nxy)
      integer nxy, isz
      integer i, ii

c  Do vector sums:
      do ii = 1, nxy, isz
        do i = ii, min(nxy,ii+isz-1)
          x(i) = y(i) 
        enddo
      enddo

      return
c
c  end of subroutine dcopy
c
      end

c************************************************************************************
c************************************************************************************
c METHOD     : matset
c DESCRIPTION: sets matrix to specified value
c ARGUMENTS  : 
c RETURNS    : 
c************************************************************************************
      subroutine matset(C, nci, ncj, val, isz)
      real*8  C(nci,ncj), val
      integer nci, ncj, isz
      integer i, ii

      do ii = 1, nci, isz
        do jj = 1, ncj, isz

          do i = ii, min(nci,ii+isz-1)
            do j = jj, min(ncj,jj+isz-1)
                C(i,j) = val
            enddo
          enddo

        enddo
      enddo

      return
c
c  end of subroutine matset
c
      end

c************************************************************************************
c************************************************************************************
c METHOD     : doop
c DESCRIPTION:
c ARGUMENTS  : 
c RETURNS    : 
c************************************************************************************
      subroutine doop(u, g, ig, nm,  nop)
      real*8  u (*), g (*)
      real*8  gtmp
      integer ig(nop*nm), nm, nop
      integer i, isz, j, k
c
c     stop "subroutine doop: do not go there!"
c  Do assignment via indirection:  
c     write(*,*), "nm==", nm
c     write(*,*), "nop=", nop
c     isz =  1
      do i = 1, nop
        gtmp = u(i)
        k    = nm*i
c     write(*,*), "g=", gtmp
c       nn   = im(i)
        do j = 1, nm
c     write(*,*), "ig=", ig(k+j)
          u(ig(k+j)+1) = gtmp
c     write(*,*), "u=", u(ig(i*nm+j))
        enddo
      enddo
c
      return
c
c  end of subroutine doop
c
      end


c************************************************************************************
c************************************************************************************
c METHOD     : isassign
c DESCRIPTION: do assignment by strided-indirection
c ARGUMENTS  : 
c RETURNS    : 
c************************************************************************************
      subroutine isassign(u, g, ig, nm, nd)
      real*8  u (*), g (*)
      real*8  gtmp
      integer nd, nm, ig(*)
      integer i, isz, j, k, m
c
c  Do assignment via indirection:  
c     write(*,*), "nm==", nm
c     write(*,*), "nop=", nop
c     isz =  1
      do i = 1, nd
        gtmp = u(i)
        k    = nm*i
c     write(*,*), "g=", gtmp
c       nn   = im(i)
        do j = 1, nm  ! stride for this i
c         write(*,*), "ig=", ig(k+j)
          m = ig(k+j)+1
          u(m) = gtmp
cc    write(*,*), "u=", u(ig(i*nm+j))
        enddo
      enddo
c
      return
c
c  end of subroutine isassign
c
      end


c************************************************************************************
c************************************************************************************
c METHOD     : isum
c DESCRIPTION: do addition using indirection 
c ARGUMENTS  : 
c RETURNS    : 
c************************************************************************************
      subroutine isum(res, u, ig, ne)
      real*8  res(*), u (*)
      integer ne, ig(*)
      integer isz, j, k
c
c  Do addition via indirection:  
c     isz =  1
      do j = 1, ne
        res(j) = 0.0;
      enddo
      do j = 1, ne
c       k = ig(j)+1
        res(j) = res(j) + u(k) 
      enddo
c
      return
c
c  end of subroutine isum
c
      end


c************************************************************************************
c************************************************************************************
c METHOD     : iprod
c DESCRIPTION: do multiplication using indirection 
c ARGUMENTS  : 
c RETURNS    : 
c************************************************************************************
      subroutine iprod(res, u, ig, ne)
      real*8  res, u (*)
      integer ne, ig(*)
      integer isz, j
c
c  Do addition via indirection:  
c     isz =  1
      res = 1.0;
      do j = 1, ne
        res = res * u(ig(j)+1)
      enddo
c
      return
c
c  end of subroutine prod
c
      end


c************************************************************************************
c************************************************************************************
c METHOD     : imax
c DESCRIPTION: find maxiumum check using indirection 
c ARGUMENTS  : 
c RETURNS    : 
c************************************************************************************
      subroutine imax(res, u, ig, ne)
      real*8  res, u (*)
      integer ne, ig(*)
      integer isz, j
c
c  Do max via indirection:  
c     isz =  1
      res = -1.0e15
      do j = 1, ne
        res = MAX( res, u(ig(j)+1) )
      enddo
c
      return
c
c  end of subroutine imax
c
      end


c************************************************************************************
c************************************************************************************
c METHOD     : imin
c DESCRIPTION: find minimum check using indirection 
c ARGUMENTS  : 
c RETURNS    : 
c************************************************************************************
      subroutine imin(res, u, ig, ne)
      real*8  res, u (*)
      integer ne, ig(*)
      integer isz, j
c
c  Do min via indirection:  
c     isz =  1
      res = 1.0e15
      do j = 1, ne
        res = MIN( res, u(ig(j)+1) )
      enddo
c
      return
c
c  end of subroutine imin
c
      end



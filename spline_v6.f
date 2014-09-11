      program main
!====================================================================
! Spline interpolation
! Comments: values of function f(x) are calculated in n base points
! then: spline coefficients are computed
!       spline interpolation is computed in 2n-1 points, 
!       a difference sum|f(u)-ispline(u)| 
!====================================================================
      implicit none
      integer, parameter :: n=8 ! base points for interpolation
      integer, parameter :: npts=1000 ! compute interpolation in npts points
      real, parameter :: H0=67.80 !Hubble parameter, Planck 2013
      real :: xmin, xmax, ymin, ymax     ! given interval of x()
      real*4,dimension (n) :: xi(n), yi(n), b(n), c(n), d(n)
      real*4,dimension(npts) :: x,xe,ys,ye,z,sumN,sumX,F,H,lines
      real*4,dimension(npts) :: xd, sumW, weight, weights
      integer i,j
      real*4 ispline, step, stepexp, total, numlin
      real*4 zstart, zend, dz
      real*4 t, k, th
      external func, gamma
      real*4 func, gamma, gauss16
!====================================================================
! open files
      open (unit=1, file='mcmc.dat')
      do i=1,n
         read (1,*,end=100) xi(i),yi(i)
      enddo
 50   format(f2.0,1x,f5.2)
 60   format(f7.4,1x,f6.2,1x,f5.2,1x,f5.2)
 100  close(1)

      xmin = 12.0
      xmax = 22.0
c      write (*,*) xmin, xmax
!====================================================================      
! step 1: generate xi and yi from  xmin, xmax, n
!      step = (xmax-xmin)/(n-1)
!      do i=1,n
!         xi(i) = xmin + step*float(i-1) 
!         yi(i) = f(xi(i)) 
!         write (*,200) xi(i), yi(i)
!      end do
!====================================================================
!  step 2: call spline to calculate spline coeficients
      call spline (xi, yi, b, c, d,n) 
!      write (*,*) d
!  step 3: interpolation at npts points
c      errav = 0.0
      step = (xmax-xmin)
      total=0.
      do i=1, npts
         xd(i)=xmin + step/(npts-1)*float(i-1)
         total=total+gauss16(gamma,0.,float(i))
      end do
      do i=1, npts
         weight(i) = 0
         sumW(i) = 0
         do j=1,i
            sumW(i) = sumW(i) + gauss16(gamma,0.0,float(j))
!                 (xd(j)**(k-1.)*exp(-xd(j)/th)/
!     +           (gauss16(gamma,0.,10000.)*th**k))
         end do
         weight(i) = 1 - sumW(i)/total
c         write (*,*) i, gauss16(gamma,0.0,float(i))
      end do
      do i=1,npts
         x(i)= xmin + step*weight(npts-i+1)
         xe(i) = 10**x(i)
         ys(i) = ispline(x(i), xi, yi, b, c, d, n)
         ye(i) = 10**ys(i)
         weights(i)=weight(npts+1-i)
! step 4: calculate quality of interpolation               
      end do
      do i=1,npts
         
         write (*,200) i, x(i), xe(i), weight(i), xe(i)!, error(i)
      end do
 200  format (i4,4d12.4)
 201  format ('           i           sum              F')    
 202  format ('           Average error',f12.5)
!====================================================================
! step 5: calculate CDF
      total=0
      do i=1,npts
         if (i.eq.npts) then 
            step=0
         else 
            step=xe(i+1)-xe(i)
         end if
         total=total+ye(i)*step
      end do
      do i=1,npts
         sumN(i)=0
         F(i)=0
         do j=1,i
            if (j.eq.i) then 
               step=0
            else 
               step=xe(j+1)-xe(j)
            endif
            sumN(i)=sumN(i)+ye(j)*step
         end do
         F(i)=sumN(i)/total
         write (*,*) i,xe(i),F(i)
      end do
!====================================================================
! step 6: calculate number of lines, n
! n = integrate [f(NHI,z) dX/dz] dNHI dz
! F = integrate [f(NHI,z) dNHI]
! dX/dz = [H0*(1+z)^2]/H(z)
! H(z) = H0 * Sqrt{[OmegaM(1+z)^3 + OmegaRAD(1+z)^4 + OmegaDE]}
!      OmegaM  = 0.27
!      OmegaDE = 0.73
!      OmegaRAD approx 0
! H = integrate [dX/dz] dz
! n = integrate [F(z) H0*(1+z)^2/H(z)] dz
!====================================================================
! calculate total absorption length
      zstart=1.87   !zstart=(wstart/1215.67)-1, wstart=3500
      zend=3.3
      dz=(zend-zstart)/npts
      total=gauss16(func,zstart,zend)
      do i=1,npts
         z(i)=zstart+dz*float(i-1)
      end do
      do i=1,npts
         sumX(i) = gauss16(func,z(1),z(i))
         H(i) = sumX(i)/total
!         write (*,*) i,sum(i),H(i)
      end do
!      do i=1,npts
!!         write (*,*) i, x(i), sumN(i), sumX(i)
!      end do
      numlin=nint(sumN(npts)*gauss16(func,zstart,zend))
      write (6,*) 'Total no. of lines = ', numlin
      write (6,*) 'NHI lower limit = ',10**xmin
      write (6,*) 'NHI upper limit = ',10**xmax


!====================================================================      
! step 6: plot data

      call PGBEGIN (0,'/xserve',2,2)
      call PGSLW(3)
      call PGENV (11.5,22.5,-26.0,-8.0,0,1)
      call PGLABEL ('log (NHI)','log f(NHI,X)','Spline in log')
      call pgsci(2)
!      call pgpt(npts,x,ys1,0)   !plot spline from MCMC
      call pgpt(n,xi,yi,-3)   !plot data from MCMC
      call pgline(npts,x,ys)
      call pgptxt(15.8,-24.0,0.0,0.5,'Prochaska spline model')

      call pgsci(1)
      call PGENV(11.5,22.5,-1e-11,2.1e-10,0,1)
      call PGLABEL ('log (NHI)','f(NHI,X)','Spline in lin')
      call pgsci(2)
      call pgline(npts,x,ye)

      call pgsci(1)
      call PGENV (11.5,22.5,0.0,1.0,0,1)!-.05,1.05,0,1)
      call PGLABEL ('log (NHI)','CDDF','Column Density Distribution
     +Function')
      call PGPTXT(14.5,0.3,0.0,0.0,'CDDF(x) =1/N \x \gS\di=1\u\ui=n\d  
     &f(x)\.\gDx')
      call PGPTXT(14.5,0.25,0.0,0.0,'N=\gS\di=1\u\ui=1000\df(x)\.\gDx')
      call pgsci(3)
      call pgline(npts,x,F)

      call pgsci(1)
      call PGENV (11.5,22.5,1.0,2500.0,0,1)
      call PGLABEL ('log NHI','X(z)','X(z)')
!      call PGLABEL ('log (NHI)','n(NHI,z) = CDDF(NHI) \x N \x X(z)',
!     +'Total number of lines')
!      call PGPTXT(14.5,800.0,0.0,0.0,'X(z) = \(2268) \dz1\u\uz2\d
!     &H\d0\u \x (1+z)\u2\d / H(z) dz')
      call pgsci(3)
      call pgline(npts,z,H)
      call pgsci(4)
!      call pgline(npts,)
      call pgend
      end program main


!======================================================================
!  Function X(z)
!======================================================================
      function func(x)
      real*4 f, x
      f = ((1+x)**2)*(0.3*((1+x)**3) + 0.7)**(-0.5)
      return
      end function func
!======================================================================
!  Function X(z)
!======================================================================
      function gamma(x)
      real*4 f, x, k
      k=1.7
      f = x**(k-1)*exp(-x)
      return
      end function gamma


!======================================================================
      subroutine spline (x, y, b, c, d, n)
!======================================================================
!  Calculate the coefficients b(i), c(i), and d(i), i=1,2,...,n
!  for cubic spline interpolation
!  s(x) = y(i) + b(i)*(x-x(i)) + c(i)*(x-x(i))**2 + d(i)*(x-x(i))**3
!  for  x(i) <= x <= x(i+1)
!  Alex G: January 2010
!----------------------------------------------------------------------
!  input..
!  x = the arrays of data abscissas (in strictly increasing order)
!  y = the arrays of data ordinates
!  n = size of the arrays xi() and yi() (n>=2)
!  output..
!  b, c, d  = arrays of spline coefficients
!  comments ...
!  spline.f90 program is based on fortran version of program spline.f
!  the accompanying function fspline can be used for interpolation
!======================================================================
      implicit none
      integer n
      real*4 x(n), y(n), b(n), c(n), d(n)
      integer i, j, gap
      real*4 h

      gap = n-1
! check input
      if ( n < 2 ) return
      if ( n < 3 ) then
         b(1) = (y(2)-y(1))/(x(2)-x(1)) ! linear interpolation
         c(1) = 0.
         d(1) = 0.
         b(2) = b(1)
         c(2) = 0.
         d(2) = 0.
         return
      end if
!
! step 1: preparation
!
      d(1) = x(2) - x(1)
      c(2) = (y(2) - y(1))/d(1)
      do i = 2, gap
         d(i) = x(i+1) - x(i)
         b(i) = 2.0*(d(i-1) + d(i))
         c(i+1) = (y(i+1) - y(i))/d(i)
         c(i) = c(i+1) - c(i)
      end do
!
! step 2: end conditions 
!
      b(1) = -d(1)
      b(n) = -d(n-1)
      c(1) = 0.0
      c(n) = 0.0
      if(n /= 3) then
         c(1) = c(3)/(x(4)-x(2)) - c(2)/(x(3)-x(1))
         c(n) = c(n-1)/(x(n)-x(n-2)) - c(n-2)/(x(n-1)-x(n-3))
         c(1) = c(1)*d(1)**2/(x(4)-x(1))
         c(n) = -c(n)*d(n-1)**2/(x(n)-x(n-3))
      end if
!
! step 3: forward elimination 
!
      do i = 2, n
         h = d(i-1)/b(i-1)
         b(i) = b(i) - h*d(i-1)
         c(i) = c(i) - h*c(i-1)
      end do
!
! step 4: back substitution
!
      c(n) = c(n)/b(n)
      do j = 1, gap
         i = n-j
         c(i) = (c(i) - d(i)*c(i+1))/b(i)
      end do
!
! step 5: compute spline coefficients
!
      b(n) = (y(n) - y(gap))/d(gap) + d(gap)*(c(gap) + 2.0*c(n))
      do i = 1, gap
         b(i) = (y(i+1) - y(i))/d(i) - d(i)*(c(i+1) + 2.0*c(i))
         d(i) = (c(i+1) - c(i))/d(i)
         c(i) = 3.*c(i)
      end do
      c(n) = 3.0*c(n)
      d(n) = d(n-1)
      end subroutine spline

!======================================================================
      function ispline(u, x, y, b, c, d, n)
!======================================================================
! function ispline evaluates the cubic spline interpolation at point z
! ispline = y(i)+b(i)*(u-x(i))+c(i)*(u-x(i))**2+d(i)*(u-x(i))**3
! where  x(i) <= u <= x(i+1)
!----------------------------------------------------------------------
! input..
! u       = the abscissa at which the spline is to be evaluated
! x, y    = the arrays of given data points
! b, c, d = arrays of spline coefficients computed by spline
! n       = the number of data points
! output:
! ispline = interpolated value at point u
!=======================================================================
      implicit none
      real*4 ispline
      integer n
      real*4  u, x(n), y(n), b(n), c(n), d(n)
      integer i, j, k
      real*4 dx

! if u is ouside the x() interval take a boundary value (left or right)
      if(u <= x(1)) then
         ispline = y(1)
         return
      end if
      if(u >= x(n)) then
         ispline = y(n)
         return
      end if

!*
!  binary search for for i, such that x(i) <= u <= x(i+1)
!*
      i = 1
      j = n+1
      do while (j > i+1)
         k = (i+j)/2
         if(u < x(k)) then
            j=k
         else
            i=k
         end if
      end do
!*
!  evaluate spline interpolation
!*
      dx = u - x(i)
      ispline = y(i) + dx*(b(i) + dx*(c(i) + dx*d(i)))
      return
      end function ispline

!==========================================================
      Function gauss16(f,a,b)
!==========================================================
! Integration of f(x) on [a,b]
! Method: Gauss 16 points  
! written by: Alex Godunov (October 2009)
! http://ww2.odu.edu/~agodunov/computing/programs/book2/Ch03/gauss.f90
!----------------------------------------------------------
! IN:
! f   - Function to integrate (supplied by a user)
! a	  - Lower limit of integration
! b	  - Upper limit of integration
! OUT:
! gauss8 - Result of integration
!==========================================================
      implicit none
      integer, parameter :: n=8
      real*4 gauss16, f, a, b
      real*4 ti(n), ci(n)
      data ti/0.0950125098, 0.2816035507, 0.4580167776, 0.6178762444,   
     &   0.7554044083, 0.8656312023, 0.9445750230, 0.9894009349/ 
      data ci/0.1894506104, 0.1826034150, 0.1691565193, 0.1495959888,
     &   0.1246289712, 0.0951585116, 0.0622535239, 0.0271524594/ 
      real*4 r, m, c
      integer i

      r = 0.0;
      m = (b-a)/2.0;
      c = (b+a)/2.0;

      do i = 1,n 
         r = r + ci(i)*(f(m*(-1.0)*ti(i) + c) + f(m*ti(i) + c))
      end do
      gauss16 = r*m
      return
      end function gauss16

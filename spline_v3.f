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
      integer, parameter :: nint=1000 ! compute interpolation in nint points
      real, parameter :: H0=67.80 !Hubble parameter, Planck 2013
      real :: xmin, xmax, ymin, ymax     ! given interval of x()
      real*4, dimension (n) :: xi(n), yi(n), b(n), c(n), d(n)
      real*4, dimension(nint) :: x,xe,ys,ye,sum,F
      integer i,j
      real*4 ispline, step, total
      real*4 z, zstart, zend, bigX
      external func, func2
      real*4 func, func2, gauss16
!====================================================================
! open files
      open (unit=1, file='mcmc.dat')
      do i=1,n
         read (1,*,end=100) xi(i),yi(i)
      enddo
 50   format(f2.0,1x,f5.2)
 60   format(f7.4,1x,f6.2,1x,f5.2,1x,f5.2)
 100  close(1)

      xmin = minval(xi)
      xmax = maxval(xi)
      write (*,*) xmin, xmax
 150  format(f5.2,2x,f6.2)
c      do i=1,n
c         write (*,150) xi(i),yi(i)
c      enddo
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
!  step 3: interpolation at nint points
c      errav = 0.0
      step = (xmax-xmin)/(nint-1)
      write(*,201) 
      do i=1, nint
         x(i) = xmin + step*float(i-1) 
         xe(i) = 10**x(i)
         ys(i) = ispline(x(i), xi, yi, b, c, d, n)
         ye(i) = 10**(ys(i))
c         write (*,200) x(i), ys(i)!, error(i)
! step 4: calculate quality of interpolation               
      end do
 200  format (3f12.4)
 201  format ('           i           sum              F')    
 202  format ('           Average error',f12.5)
!====================================================================
! step 5: calculate CDF
      total=0
      do i=1,nint
         total=total+ye(i)
      end do
      do i=1,nint
         sum(i)=0
         F(i)=0
         do j=1,i
            sum(i)=sum(i)+ye(j)
         end do
         F(i)=sum(i)/total
         write (*,*) i,sum(i),F(i)
      end do
!====================================================================
! step 6: calculate number of lines, n
! n = integrate [f(NHI,z) dX/dz] dNIH dz
! F = integrate [f(NHI,z) dNHI]
! dX/dz = [H0*(1+z)^2]/H(z)
! H(z) = H0 * Sqrt{[OmegaM(1+z)^3 + OmegaRAD(1+z)^4 + OmegaDE]}
!      OmegaM  = 0.27
!      OmegaDE = 0.73
!      OmegaRAD approx 0
! H = integrate [dX/dz] dz
! n = integrate [F(z) H0*(1+z)^2/H(z)] dz
      zstart=1.700
      zend=2.141
      bigX = 0
      bigX = gauss16(func,zstart,zend)
      write (*,*) 'X of func from gauss16:',bigX
      bigX = gauss16(func2,zstart,zend)
      write (*,*) 'X of func2 from gauss16:',bigX
!====================================================================      
! step 6: plot data

      call PGBEGIN (0,'/xserve',2,2)
      call PGENV (11.5,22.5,-26.0,-8.0,0,1)
      call PGLABEL ('log (NHI)','log f(NHI,X)','Spline in log')
      call PGSLW(3)
      call pgsci(2)
!      call pgpt(nint,x,ys1,0)   !plot spline from MCMC
      call pgpt(n,xi,yi,-3)   !plot data from MCMC
      call pgline(nint,x,ys)
      call pgptxt(15.8,-24.0,0.0,0.5,'Prochaska spline model')

      call pgsci(1)
      call PGENV(11.5,22.5,-1e-11,2.1e-10,0,1)
      call PGLABEL ('log (NHI)','f(NHI,X)','Spline in lin')
      call pgsci(2)
      call pgline(nint,x,ye)

      call pgsci(1)
      call PGENV (11.5,22.5,-.05,1.05,0,1)
      call PGLABEL ('log (NHI)','CDF','CDF')
      call pgsci(3)
      call pgline(nint,x,F)

      call pgsci(1)
      call PGENV (1e11,1e22,-1e-11,2.1e-10,0,1)
      call PGLABEL ('NHI','f(NHI)','Test')
      call pgsci(2)
      call pgline(nint,xe,ye)
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
      function func2(x)
      real*4 f, x, x1,x2
      x1 = 0.3*(1+x)**2*(x+1/1.35)
      x2 = 0.7*((1+x)**2/0.68 + x**2 + 2*x)
      f = (1+x)**2*(x1+x2)**(-0.5)
      return
      end function func2


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

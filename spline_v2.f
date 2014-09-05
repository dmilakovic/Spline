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
      real*4, dimension(nint) :: x,y,ys,ye,sum,F
      integer i,j
      real*4 ispline, step, total, z, zend, H

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
      do i=1,n
         write (*,150) xi(i),yi(i)
      enddo
      
! step 1: generate xi and yi from  xmin, xmax, n
!      step = (xmax-xmin)/(n-1)
!      do i=1,n
!         xi(i) = xmin + step*float(i-1) 
!         yi(i) = f(xi(i)) 
!         write (*,200) xi(i), yi(i)
!      end do

!  step 2: call spline to calculate spline coeficients
      call spline (xi, yi, b, c, d,n) 
!  step 3: interpolation at nint points
c      errav = 0.0
      step = (xmax-xmin)/(nint-1)
      write(*,201) 
      do i=1, nint
         x(i) = xmin + step*float(i-1) 
         y(i) = f(x) 
         ys(i) = ispline(x(i), xi, yi, b, c, d, n)
         ye(i) = 10**(ys(i))
c         write (*,200) x(i), ys(i)!, error(i)
! step 4: calculate quality of interpolation               
      end do
 200  format (3f12.4)
 201  format ('        x        spline1      splin2')    
 202  format ('           Average error',f12.5)
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
      zend=3
      call qgaus(0,zend,f,H)

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
      call pgend
      end program main


!======================================================================
!  Function f(z)
!======================================================================
      function f(z)
      double precision f, z
      f = (1+z)^2/(0.27*(1+z)^3 + 0.73)^(0.5)
      end function f

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
      end function ispline

!=======================================================================
      SUBROUTINE qgaus(func,a,b,ss)
!======================================================================
! INPUT: 
! a, b = integration limits
! func = function to be integrated
! OUTPUT:
! ss = integrand
!=======================================================================
      IMPLICIT NONE
      REAL a,b,ss,func
      EXTERNAL func
!          Returns as ss the integral of the function func between a and b, by              ten-point Gauss-Legendre integration: the function is evaluated                 exactly ten times at interior points in the range of integration.
      INTEGER j

      REAL dz,zm,zr,w(5),zx(5) !The abscissas and weights.
      SAVE w, zx
      DATA w/.2955242247,.2692667193,.2190863625,.1494513491,
     +      .0666713443/
      DATA zx/.1488743389,.4333953941,.6794095682,.8650633666,
     +      .9739065285/
      zm=0.5*(b+a)
      zr=0.5*(b-a)
      ss=0 !Will be twice the average value of the function, since the ten
!           weights (five numbers above each used twice) sum to 2.
      do j=1,5 
         dz=zr*zx(j)
         ss=ss+w(j)*(func(zm+dz)+func(zm-dz))
      enddo
      ss=zr*ss !Scale the answer to the range of integration.
      return
      END

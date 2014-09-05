      program main
!====================================================================
! Spline interpolation
! Comments: values of function f(x) are calculated in n base points
! then: spline coefficients are computed
!       spline interpolation is computed in 2n-1 points, 
!       a difference sum|f(u)-ispline(u)| 
!====================================================================
      implicit none
      integer, parameter :: n1=8 ! base points for interpolation
      integer, parameter :: n2=37
      integer, parameter :: nint=1000 ! compute interpolation in nint points
      real :: xmin, xmax, ymin, ymax     ! given interval of x()
      real*4, dimension (n1) :: xi1(n1), yi1(n1), b1(n1), c1(n1), d1(n1)
      real*4, dimension (n2) :: xi2(n2), yi2(n2), b2(n2), c2(n2), d2(n2)
      real*4, dimension (n2) :: errup2(n2), errlow2(n2)
      real*4, dimension(nint) :: x(nint),ys1(nint), ys2(nint)
      real :: z(nint), zs(nint)
      integer i
      real*4 f, ispline, step

! open files
      open (unit=1, file='mcmc.dat')
      do i=1,n1
         read (1,*,end=100) xi1(i),yi1(i)
      enddo
      open (unit=2, file='data.dat')
      do i=1,n2
         read (2,*,end=110) xi2(i),yi2(i),errup2(i),errlow2(i)
      enddo
 50   format(f2.0,1x,f5.2)
 60   format(f7.4,1x,f6.2,1x,f5.2,1x,f5.2)
 100  close(1)
 110  close(2)

      xmin = 10.0
      xmax = 23.0
      write (*,*) '----- file 1 -----'
 150  format(f5.2,2x,f6.2)
 160  format(f7.4,2x,f6.2,2x,f5.2,2x,f5.2)
      do i=1,n1
         write (*,150) xi1(i),yi1(i)
      enddo
      write (*,*) '----- file 2 -----'
      do i=1,n2
         write (*,160) xi2(i),yi2(i),errup2(i),errlow2(i)
      enddo
! step 1: generate xi and yi from f(x), xmin, xmax, n
!      step = (xmax-xmin)/(n-1)
!      do i=1,n
!         xi(i) = xmin + step*float(i-1) 
!         yi(i) = f(xi(i)) 
!         write (*,200) xi(i), yi(i)
!      end do

!  step 2: call spline to calculate spline coeficients
      call spline (xi1, yi1, b1, c1, d1,n1) 
      call spline (xi2, yi2, b2, c2, d2,n2)
!  step 3: interpolation at nint points
c      errav = 0.0
      step = (xmax-xmin)/(nint-1)
      write(*,201) 
      do i=1, nint
         x(i) = xmin + step*float(i-1) 
c         y(i) = f(x) 
         ys1(i) = ispline(x(i), xi1, yi1, b1, c1, d1, n1)
         ys2(i) = ispline(x(i), xi2, yi2, b2, c2, d2, n2)
c         error(i) = ys(i)-y
         write (*,200) x(i), ys1(i), ys2(i)!, error(i)
! step 4: calculate quality of interpolation               
      end do
 200  format (3f12.4)
 201  format ('        x        spline1      splin2')    
 202  format ('           Average error',f12.5)
! step 5: plot data
      ymin=0.0
      ymax=0.0
      do i=1,nint
        if(ys1(i).gt.ymax)ymax=ys1(i)
      end do
      ymin=ymax
      do i=1,nint
        if(ys1(i).lt.ymin)ymin=ys1(i)
      end do
      ymax=ymax-0.1*ymax
      ymin=ymin+0.1*ymin
      xmax=1.1*maxval(x)
      xmin=0.9*minval(x)

      call PGBEGIN (0,'/vcps',1,1)
      call PGENV (11.5,22.5,-26.0,-8.0,1,1)
      call PGLABEL ('log (NHI)','log f(NHI,X)','Spline')
      call PGSLW(3)
      call pgsci(2)
!      call pgpt(nint,x,ys1,0)   !plot spline from MCMC
      call pgpt(n1,xi1,yi1,-3)   !plot data from MCMC
      call pgline(nint,x,ys1)
      call pgptxt(15.8,-24.0,0.0,0.5,'Prochaska spline model')
      call pgsci(3)
!      call pgpt(nint,x,ys2,0)   !plot spline from whole sample
      call pgpt(n2,xi2,yi2,-4)   !plot data from whole sample
!      call pgline(nint,x,ys2)
      call pgptxt(15.8,-23.5,0.0,0.5,'K13R13 & N12 & OPB07')
      call pgend
      end program main


!
!  Function f(x)
!
      function f(x)
      double precision f, x
      f = sin(x) 
      end function f

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

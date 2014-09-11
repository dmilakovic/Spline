#makefile for spline
FC=gfortran

pgplot=-L/opt/local/lib -lpgplot
X11=-L/usr/X11/lib -lX11
input=spline_v6.f
output=spline_v6

all:
	$(FC) -o $(output) $(input) $(pgplot) $(X11)

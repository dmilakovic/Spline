#makefile for spline
FC=gfortran

pgplot=-L/opt/local/lib -lpgplot
X11=-L/usr/X11/lib -lX11
input=spline_v4.f
output=spline_v4

all:
	$(FC) -o $(output) $(input) $(pgplot) $(X11)

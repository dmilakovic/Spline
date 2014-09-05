import numpy as np
import matplotlib.pyplot as plt
from scipy import interpolate
import random as rnd

f1=open('mcmc.dat','r')
i=0
x1=np.empty(8)
y1=np.empty(8)
for line in f1:
    values1=line.split(',')
    x1[i]=float(values1[0])
    y1[i]=float(values1[1])
    i=i+1
f1.close()

f2=open('data.dat','r')
i=0
x2=np.empty(37)
y2=np.empty(37)
for line in f2:
    values2=line.split(',')
    x2[i]=float(values2[0])
    y2[i]=float(values2[1])
    i=i+1
f2.close()

x=np.append(x1,x2)
xmin=min(x)
xmax=max(x)

npoints=1000
dx=(xmax-xmin)/npoints
tck = interpolate.splrep(x1,y1,s=0)
xnew = np.arange(xmin,xmax,dx)
ynew = interpolate.splev(xnew,tck,der=0)

S=sum(ynew)
F=np.empty(npoints)
for i in xrange(1,npoints):    
    F[i]=sum(ynew[0:i])/S

plt.subplot(2,1,1)
plt.plot(xnew, ynew, '-r', linewidth=2)
plt.scatter(x1, y1, c='r', marker='v')
plt.scatter(x2, y2, c='g', marker='o')
plt.legend(['Spline as per P13', 'P13 parameters', 'P13 data'])
plt.ylabel('log f(N,X)')
plt.xlim(11.5, 22.5)
plt.ylim(-27, -9)
plt.title('Cubic-spline interpolation')

plt.subplot(2,1,2)
plt.plot(xnew,F)
plt.xlim(11.5,22.5)
plt.ylim(-.05,1.05)
plt.ylabel('F(N)')
plt.xlabel('log N')
plt.show()



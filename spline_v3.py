import numpy as np
import matplotlib.pyplot as plt
from scipy import interpolate, integrate
import random as rn

#-------------------------------------------------------
#                read data from files
f1=open('mcmc.dat','r')
i=0
xMCMC=np.empty(8,dtype=np.float64)
yMCMC=np.empty(8,dtype=np.float64)
yexp=np.empty(8,dtype=np.float64)
for line in f1:
    values1=line.split(',')
    xMCMC[i]=float(values1[0])
    yMCMC[i]=float(values1[1])
    i=i+1
f1.close()

f2=open('data.dat','r')
i=0
xDATA=np.empty(37,dtype=np.float64)
yDATA=np.empty(37,dtype=np.float64)
for line in f2:
    values2=line.split(',')
    xDATA[i]=float(values2[0])
    yDATA[i]=float(values2[1])
    i=i+1
f2.close()

x=np.append(xMCMC,xDATA)
xmin=min(x)
xmax=max(x)
#-------------------------------------------------------
#                  interpolate spline
npoints=1000
dx=(xmax-xmin)/npoints
tck1 = interpolate.splrep(xMCMC,yMCMC,s=0)
NHI = np.arange(xmin,xmax,dx)
logf = interpolate.splev(NHI,tck1,der=0)
f = 10**logf
#-------------------------------------------------------
#   determine Cumulative Distribution Function (CFD)
Sumf=float(sum(f))
count=0
CDF=np.empty(npoints,dtype=np.float64)
for i in xrange(0,npoints):    
    CDF[i]=sum(f[0:i])/Sumf
    if CDF[i] == 1.0:
        count=count+1
#    print i,'%.25e'%CDF[i]
    print i,'%.5e'%sum(f[0:i]),'%.25e'%CDF[i]
print 'times that value 1.0 is repeated :', count

#interpolate CDF
tck2 = interpolate.splrep(NHI,CDF,s=0)
#-------------------------------------------------------
#  determine X = integrate (dX/dz) dz
def integrand(x):
    return (1+x)**2/np.sqrt(0.3*(1+x)**3 + 0.7)
X = integrate.quad(integrand,1.9,3.2)
print 'X:',X
#-------------------------------------------------------
#                  determine mean
epsilon=0.005
dummy,count=0,0
for i in xrange(0,npoints):
    if (abs(CDF[i]-0.5) <= epsilon) and (abs(CDF[i]-0.5) < abs(CDF[dummy]-0.5)):
        dummy=i
        count=count+1

xmean1=NHI[dummy]
ymean1=CDF[dummy]
print 'number of elements within |CDF[i]-0.5|<=',epsilon,':',count
sigma1=np.sqrt(sum((NHI-xmean1)**2)/npoints)
xmean2=sum(NHI*f)/Sumf
mom2=sum(NHI**2*f)/Sumf
sigma2=np.sqrt(mom2-xmean2**2)
ymean2=interpolate.splev(xmean2,tck2,der=0)
print 'mean1 from |F[i]-0.5| <=',epsilon,':', xmean1
print 'mean2 from E(x)=1/S sum(x*f(x)):', xmean2
print 'F(mean1), read from matrix:', ymean1
print 'F(mean2), interpolated:', ymean2

#-------------------------------------------------------    
#                    plot data
fig1=plt.figure()
sub1=fig1.add_subplot(2,1,1)
sub1.plot(NHI, logf, '-r', linewidth=2)
sub1.scatter(xMCMC, yMCMC, c='r', marker='v')
sub1.scatter(xDATA, yDATA, c='g', marker='o')
sub1.legend(['Spline as per P13', 'P13 parameters', 'P13 data'])
sub1.set_ylabel('log [f(N,X)]')
sub1.set_xlim(11.5, 22.5)
sub1.set_ylim(-27, -9)
#axis1.yscale('linear')
#fig1.title('Cubic-spline interpolation')

sub2=fig1.add_subplot(2,1,2)
sub2.plot(NHI,CDF)
sub2.set_xlim(11.5,22.5)
sub2.set_ylim(-.05,1.05)
sub2.set_ylabel('CDF')
sub2.set_xlabel('log N')
plt.savefig('fig1.pdf')

fig2=plt.figure()
sub3=fig2.add_subplot(2,1,1)
sub3.plot(NHI, f, '-r', linewidth=2)
sub3.set_ylabel('f(N)')
sub3.set_xlim(11.5,22.5)
sub3.set_ylim(-1e-11,2.1e-10)

sub4=fig2.add_subplot(2,1,2)
sub4.plot(NHI,CDF)
sub4.set_xlim(11.5,22.5)
sub4.set_ylim(-.05,1.05)
sub4.set_ylabel('CDF')
sub4.set_xlabel('log N')
plt.savefig('fig2.pdf')


fig3=plt.figure()
ax=fig3.add_subplot(111)
ax.set_title('Histogram of sampled data')
n,bins,patches=ax.hist(NHI,30,normed=1,facecolor='green',alpha=0.75)
ax.set_xlabel('NHI')
ax.set_ylabel('Number')
plt.savefig('fig3.pdf')

## Field free molecular alignment
# Schrodinger Equation Implementation
# Tragically slow
# Craig Benko, 2014.07.31

# General libraries
# from pylab import *
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import ode
import time
import math
import seaborn as sns
from phi_st import phi_st
#import SI2au as SI2au

#配色
sns.set();





# My libraries, if using ipython, restart if these are modified.
from boltzmann import Boltzmann
from EA import EA
from cosfun import *
from molecules import *
from const import *
from SI2au import SI2au

# close old plots.
plt.close('all')
plt.style.use('default');

# start time
timer = time.time()

## tuneable parameters
molecule = 'N2'
pulse_FWHM = 100e-15 #FWHM duration of the intesity of the pulse in seconds
# I = .1 #in 10**14 W cm**-2
TemperatureK = 20  #in Kelvin


##Calculated Parameters
#molecular parameters
B = B_dict[molecule]*1.98648e-23 #rotational constant in ground state in Joules
D = D_dict[molecule]*1.98648e-23 #centrifugal distorion in ground state in Joules
delta_alpha =  d_alpha_dict[molecule] #anisotropic polarizability
Jmax = Jmax_dict[molecule] #approximate max J
Temperature = k*TemperatureK/B
Jweight = Boltzmann(Temperature, 70, molecule) #Boltzmann distribution

#laser parameters
# old electric field amplitude
# E0 = 2.74*10**10*I**.5 # old electric field amplitude

sigma = pulse_FWHM*B/hbar
I = 5e18; # 单位：W / m^2 = kg / s^3  # 光强
I = SI2au(I,'I'); # a.u.
E0 = math.sqrt(2*I/(epsilon0*c)); # electric field amplitude
strength=0.5*4*np.pi*epsilon0*delta_alpha*E0**2/B #in rotational constan
lambda_m = 400e-9  # 波长，单位：米# 定义 lambda
lambda_au = SI2au(lambda_m, 'a0');# 调用 SI2au 函数进行单位转换
c = 137.0359990838876; # a.u. # 光速
epsilon0 = 0.079577471545795; # a.u.# 真空中介电常数
T = lambda_au/c #光周期
TOAM = 2 # 拓扑荷
nu0 = 1./(lambda_au/c);   # a.u.# 时间频率
omega = nu0 * 2*math.pi; # a.u.# 时间角频率



# Cartesian coordinate system and time


#精度设置

numxx=300;
xxr=12
xx = np.linspace(-1*xxr, xxr, numxx);  
# 从-12到12，总共有(numx)个点

# x = np.arange(-12, 12, 0.03)  从-12开始，到12结束（包含12），步长为0.03
tl = 8;#时域长度(tl为脉冲周期个数，实际的时间长度为tl*T周期)
numt=200;#时域采样点个数
delta = 2*tl/numt  #时间间隔 delta*T 最好小于等于 0.1 a.u.
tt =  np.arange(0, 2*tl*T,delta*T) 
dt = delta*T

# 输出参数
print("波长（原子单位）:", lambda_au,"电场强度",E0)




## RHS of the Schrodinger Equation
def rhs(t,x,Mparm,xxs):
    dx = np.array(zeros(Jmax, dtype = 'complex'))
    # Delta_omega = pulse(t, strength, sigma)
    Delta_omegatotal = EA(xxs, t, c, lambda_au, omega, TOAM)
    Delta_omega=Delta_omegatotal[0];

    for k in range(Jmax):
        if k == 0 or k == 1:
            dx[k] = -1j*(x[k]*(k*(k+1) - D/B*k**2*(k+1)**2 - Delta_omega) -
                Delta_omega*x[k]*c2(k,Mparm) -
                Delta_omega*x[k+2]*cp2(k,Mparm))
        elif k == Jmax - 2 or k == Jmax-1:
            dx[k] = -1j*(x[k]*(k*(k+1) - D/B*k**2*(k+1)**2 - Delta_omega) -
                Delta_omega*x[k-2]*cm2(k,Mparm) -
                Delta_omega*x[k]*c2(k,Mparm))
        else:
            dx[k] = -1j*(x[k]*(k*(k+1) - D/B*k**2*(k+1)**2 - Delta_omega) -
                Delta_omega*x[k-2]*cm2(k,Mparm) -
                Delta_omega*x[k+2]*cp2(k,Mparm) -
                Delta_omega*x[k]*c2(k,Mparm))
    return dx

## Initialize
tend = 2*sigma; #dt = .04*sigma
#tt = np.linspace(0,5,1000)
# tt = np.arange(-M, M, delta) 
cos2 = np.zeros(np.shape((tt.size,xx.size)),dtype = 'complex')
Cstor = np.zeros((Jmax,int(2*Jmax+1), Jmax,numxx), dtype = 'complex')
start = np.zeros(Jmax, dtype = 'complex')




qx=np.zeros(np.shape(xx))
qt=np.zeros(np.shape(tt))
cx=np.zeros((np.size(xx),np.size(tt)))
phite=np.zeros((np.size(xx),np.size(tt)))
tphi=tt-tt[-1]/2

X3,Y3=np.meshgrid(tphi,xx)

for y in range(numt):
    phite[:,y]=phi_st(xx,tt[y],tt)
    cx[:,y]=EA(xx, tt[y], c, lambda_au, omega, TOAM,tt)[0]
    
fig1=plt.figure()
ax1=fig1.add_subplot(221)
ax2=fig1.add_subplot(222)
ax3=fig1.add_subplot(223)
ax4=fig1.add_subplot(224)
ax1.contourf(X3, Y3, cx, vmin=0)

ax2.plot(tt,cx[20,:])
ax3.plot(xx,cx[:,20])


h1=ax4.contourf(X3, Y3, phite, vmin=0)
cb=plt.colorbar(h1)





import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
sns.set()

def cartesian_to_polar(x, y):
    r = np.sqrt(x**2 + y**2)
    theta = np.arctan2(y, x)
    return r, theta






N = 200;
lambdaz = 632e-9; #%波长为632nm
k = 2*np.pi/lambdaz; #%波数
w0 = 3; #%束腰半径
x = np.linspace(-100,100,N);
y = np.linspace(-100,100,N);
X,Y= np.meshgrid(x,y);
r,theta=cartesian_to_polar(X,Y)
beta = 50*np.pi/180;
fig1=plt.figure();
l =6;
E1 = (r/w0)**abs(l)*np.exp(-r**2/w0**2)*np.exp(-1j*l*theta);
ax1=fig1.add_subplot(221)


h1=ax1.contourf(X,Y,np.abs(E1)**2,200);
cb1=plt.colorbar(h1);
ax2=fig1.add_subplot(222)
h2=ax2.contourf(X,Y,np.angle(E1));
cb2=plt.colorbar(h2);
ax3=fig1.add_subplot(223,projection='3d')

ax3.plot_surface(X,Y,np.abs(E1)**2) #%三维





from stov import stov
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
sns.set()
xx=np.linspace(-200,200,400)
tt=np.linspace(-200,200,400)
qx=np.zeros(np.shape(xx))
qt=np.zeros(np.shape(tt))
cx=np.zeros((np.size(xx),np.size(tt)))
ampx=np.zeros((np.size(xx),np.size(tt)))
px=np.zeros((np.size(xx),np.size(tt)))
pstx=np.zeros((np.size(xx),np.size(tt)))
for j in range(np.size(xx)):
    for y in range(np.size(tt)):
        cx[j,y]=np.real(stov(tt[y],xx[j],TOAM=1,ttm=300)[0])
        ampx[j,y]=stov(tt[y],xx[j],TOAM=1,ttm=300)[1]
        px[j,y]=stov(tt[y],xx[j],TOAM=1,ttm=300)[2]
        pstx[j,y]=stov(tt[y],xx[j],TOAM=1,ttm=300)[3]
        
X,Y= np.meshgrid(xx,tt);
fig1=plt.figure()
ax1=fig1.add_subplot(221)
h1=ax1.contourf(X,Y,cx,40);
ax2=fig1.add_subplot(222)
h2=ax2.contourf(X,Y,px,200);
ax3=fig1.add_subplot(223)
h3=ax3.contourf(X,Y,ampx,200);
ax4=fig1.add_subplot(224)
h4=ax4.contourf(X,Y,pstx,200,cmap='RdBu')





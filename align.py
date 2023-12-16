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
#import SI2au as SI2au

#配色
sns.set();




# My libraries, if using ipython, restart if these are modified.
from boltzmann import Boltzmann
from cosfun import *
from molecules import *
from const import *
from SI2au import SI2au
from stov import stov



# close old plots.
plt.close('all')
plt.style.use('default');

# start time
timer = time.time()

## tuneable parameters
molecule = 'N2'


#pulse_FWHM = 200e-15 #FWHM duration of the intesity of the pulse in seconds 单位：s
# I=0.5e18  # 单位：W / m^2 = kg / s^3
I = .1#in 10**14 W cm**-2
TemperatureK = 80  #in Kelvin
E0 = 2.74*10**10*I**.5 # electric field amplitude
lambda_m = 400e-9  # 波长，单位：米# 定义 lambda

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



# I = 5e18; # 单位：W / m^2 = kg / s^3  # 光强
# I = SI2au(I,'I'); # a.u.
# E0 = np.sqrt(2*I/(epsilon0*c)); # electric field amplitude

# I1 = SI2au(I,'I');#光强 #a.u.
# E0=np.sqrt(2*I1/(epsilon0*c))#a.u.



# lambda_au = SI2au(lambda_m, 'a0');# 调用 SI2au 函数进行单位转换
c=2.9979248e8;
#c = 137.0359990838876; # a.u. # 光速
epsilon0 = 8.854187817e-12;
#epsilon0 = 0.079577471545795; # a.u.# 真空中介电常数
strength=2*4*np.pi*epsilon0*delta_alpha*E0**2/B #in rotational constant
# strength=0.5*4*np.pi*epsilon0*delta_alpha*E0**2/B #in rotational constant


T = lambda_m/c #光周期 单位s
TOAM = 2 # 拓扑荷
nu0 = 1/T;   # 时间频率
omega = nu0 * 2*np.pi; # 时间角频率


pulse_FWHM=150*T;
sigma = pulse_FWHM*B/hbar;


# Cartesian coordinate system and time


#精度设置

numxx=15;
xxr=2.4e-4;
xx = np.linspace(-1*xxr, xxr, numxx);  
# 从-12到12，总共有(numx)个点

# x = np.arange(-12, 12, 0.03)  从-12开始，到12结束（包含12），步长为0.03
#tl = 12;#时域长度(tl为脉冲周期个数，实际的时间长度为tl*T周期)
#numt=300;#时域采样点个数
#delta = 2*tl/numt  #时间间隔 delta*T 最好小于等于 0.1 a.u.
#tt =  np.arange(-2*tl*T, 2*tl*T,delta*T) 
#dt = delta*T

# 输出参数
#print("波长（原子单位）:", lambda_au,"电场强度",E0)
print("波长:", lambda_m,"电场强度",E0)

## Initialize
tend = 3*sigma; 
dt = 0.04*sigma;
tt=np.linspace(0,6,1000);
ttmax=0

# stovp = stov(t,y,strength,TOAM,lambdaz,ttm)





## RHS of the Schrodinger Equation
def rhs(t,x,y,Mparm,strength,TOAM,lambdaz,ttm,Bm):
    dx = np.array(zeros(Jmax, dtype = 'complex'))
    # Delta_omega = pulse(t, strength, sigma)
    Delta_omegatotal = stov(t,y,strength,TOAM,lambdaz,ttm,Bm)
    Delta_omega=Delta_omegatotal[1];

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




# tt = np.arange(-M, M, delta) 
cos2 = np.zeros((tt.size,xx.size),dtype = 'complex')
Cstor = np.zeros((Jmax,int(2*Jmax+1), Jmax,numxx), dtype = 'complex')
start = np.zeros(Jmax, dtype = 'complex')





## Integrate Schrodinger Eq. Loop over all initial wavefunctions |J,M>
for R in range(numxx):
    for J in range(Jmax):
        for M in range(J+1):
            #create ODE
            s = ode(rhs).set_f_params(xx[R],M,strength,TOAM,lambda_m,ttmax,B).set_integrator('zvode',atol = 1e-5,rtol = 1e-4, order = 9)
            #initialize
            init = 0*start
            init[J] = 1
            s.set_initial_value(init.tolist(),0)
            solnt = []
            solny = []
            #integrate
            while s.successful() and (s.t < tend):
                sst=s.t+dt
                s.integrate(sst);
                solnt.append(s.t)
                solny.append(s.y/np.sum(s.y*np.conj(s.y))**.5)
        #store
            Cstor[J,M,:,R] = np.transpose(solny)[:,-1] 
            Cstor[J,-1*M,:,R] = np.transpose(solny)[:,-1]

## Expectation value, incoherent, thermal average.
for R in range(numxx):
    for J in range(Jmax):
        for M in range(-J,J+1):
            for jj in range(Jmax-2):
                w = 4*jj+6
                phi = np.angle(Cstor[J,M,jj,R])-np.angle(Cstor[J,M,jj+2,R])
                cos2[:,R]+= Jweight[J]/(2*J+1)*(abs(Cstor[J,M,jj,R])**2*c2(jj,M)+abs(Cstor[J,M,jj,R])*abs(Cstor[J,M,jj+2,R])*cp2(jj,M)*np.cos(w*tt+phi))

## End program
elapsed = time.time() - timer
print('\n Program took ' + str(round(elapsed)) + ' s to run. \n')
print('\n' + molecule+' at '+ str(I) + ' x 10$^{14}$ W cm$^{-2}$ at ' + str(TemperatureK) + ' K\n')

#Plot result, <cos**2\theta>
import math
plt.close('all');
fig1=plt.figure();

samp=7;#取样点个数
pa=80
new_blues=sns.color_palette("Blues_r",int(pa))[0:int(pa):int(pa/samp)+1];
sns.set_palette(new_blues);
sns.color_palette(new_blues)


ax1=fig1.add_subplot(121)
ax2=fig1.add_subplot(222)
ax7=fig1.add_subplot(224)
Y,Z=np.meshgrid(tt*hbar/B*10**12,xx)








ax2.contourf(Y,Z,np.transpose(np.real(cos2)),60)

for plotnum in range(0,numxx,math.ceil(numxx/samp)):
    ax1.plot(tt*hbar/B*10**12,np.real(cos2[:,plotnum])+np.max(cos2*(plotnum/numxx)))
    ax2.axhline(xx[plotnum],color='w',linestyle='--',linewidth=1.5)




cx1=np.zeros((np.size(xx),np.size(tt)))
ampx1=np.zeros((np.size(xx),np.size(tt)))
px1=np.zeros((np.size(xx),np.size(tt)))
pstx1=np.zeros((np.size(xx),np.size(tt)))




numsx=400
xxz=np.linspace(-xxr,xxr,numsx)
ttz=np.linspace(-3*sigma,3*sigma,400)
ttm=np.max(ttz)

for j in range(np.size(xx)):
    for y in range(np.size(ttz)):
        cx1[j,y]=np.real(stov(ttz[y],xx[j],TOAM=1,ttm=ttmax)[0])
        ampx1[j,y]=stov(ttz[y],xx[j],TOAM=1,ttm=ttmax)[1]
        px1[j,y]=stov(ttz[y],xx[j],TOAM=1,ttm=ttmax)[2]
        pstx1[j,y]=stov(ttz[y],xx[j],TOAM=1,ttm=ttmax)[3]

ax7.contourf(Y,Z,ampx1,200)
ax7.set_xlim(0,7)

for plotnum in range(0,numxx,math.ceil(numxx/samp)):
    ax7.axhline(xx[plotnum],color='w',linestyle='--',linewidth=1.5)
# ax7.set_xlim(0,0.1)





numsx=400
xxz=np.linspace(-xxr,xxr,numsx)
ttz=np.linspace(-3*sigma,3*sigma,400)
ttm=np.max(ttz)

qx=np.zeros(np.shape(xxz))
qt=np.zeros(np.shape(ttz))
cx=np.zeros((np.size(xxz),np.size(ttz)))
ampx=np.zeros((np.size(xxz),np.size(ttz)))
px=np.zeros((np.size(xxz),np.size(ttz)))
pstx=np.zeros((np.size(xxz),np.size(ttz)))
for j in range(np.size(xxz)):
    for y in range(np.size(ttz)):
        cx[j,y]=np.real(stov(ttz[y],xxz[j],TOAM=1,ttm=0)[0])
        ampx[j,y]=stov(ttz[y],xxz[j],TOAM=1,ttm=0)[1]
        px[j,y]=stov(ttz[y],xxz[j],TOAM=1,ttm=0)[2]
        pstx[j,y]=stov(ttz[y],xxz[j],TOAM=1,ttm=0)[3]
        
U,V= np.meshgrid(xxz,ttz+ttm);
fig2=plt.figure()
ax3=fig2.add_subplot(221)
h1=ax3.contourf(V,U,cx,40);
ax4=fig2.add_subplot(222)
h2=ax4.contourf(V,U,px,200);
ax5=fig2.add_subplot(223)
h3=ax5.contourf(V,U,ampx,200);
ax6=fig2.add_subplot(224)
for hi in range(numxx):
    h4=ax6.plot(ttz+ttm,ampx[int(numsx/numxx*hi)])
    #ax3.axhline(xxz[int(numsx/nmuxx*hi)],color='w',linestyle='--',linewidth=1.5)
    #ax5.axhline(xxz[int(numsx/numxx*hi)],color='w',linestyle='--',linewidth=1.5)
h4=ax6.plot(ttz+ttm,ampx[int(numsx/2)])



#plt.plot(tt*hbar/B*10**12,np.real(cos2),'k-')
#plt.plot()
# plt.xlabel('Time [ps]')
# plt.ylabel('<cos$^2\Theta$>')
#plt.title(molecule+' at '+ str(I) + ' x 10$^{14}$ W cm$^{-2}$ at ' + str(TemperatureK) + ' K')
#plt.grid()
#plt.ylim(0,1)
#plt.show()

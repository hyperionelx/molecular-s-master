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
#import SI2au as SI2au

# My libraries, if using ipython, restart if these are modified.
from boltzmann import Boltzmann
from EA import EA
from cosfun import *
from molecules import *
from const import *
from SI2au import SI2au

# close old plots.
plt.close('all')

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
TOAM = 1 # 拓扑荷
nu0 = 1./(lambda_au/c);   # a.u.# 时间频率
omega = nu0 * 2*math.pi; # a.u.# 时间角频率



# Cartesian coordinate system and time
# x = np.linspace(-12, 12, num=int((12 - (-12)) / 0.03) + 1)  # 从-12到12，总共有(num)个点
x = np.arange(-12, 12, 0.03)  # 从-12开始，到12结束（包含12），步长为0.03
M = 12
delta = 2*M/3000  #时间间隔 delta*T 最好小于等于 0.1 a.u.
dt = delta*T
t =  np.arange(-M, M, delta) 
t = t*T





# 输出参数
print("波长（原子单位）:", lambda_au,"电场强度",E0)




## RHS of the Schrodinger Equation
def rhs(t, x, Mparm):
    dx = np.array(zeros(Jmax, dtype = 'complex'))
    # Delta_omega = pulse(t, strength, sigma)
    Delta_omega = EA(x, t, c, lambda_au, omega, TOAM)

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
dt = delta*T
#tt = np.linspace(0,5,1000)
tt = np.arange(-M, M, delta) 
cos2 = np.zeros(tt.size,dtype = 'complex')
Cstor = np.zeros((Jmax,int(2*Jmax+1), Jmax), dtype = 'complex')
start = np.zeros(Jmax, dtype = 'complex')

## Integrate Schrodinger Eq. Loop over all initial wavefunctions |J,M>
for J in range(Jmax):
    for M in range(J+1):
        #create ODE
        s = ode(rhs).set_f_params(M).set_integrator('zvode',atol = 1e-5,rtol = 1e-4, order = 9)
        #initialize
        init = 0*start
        init[J] = 1
        s.set_initial_value(init.tolist(),0)
        solnt = []
        solny = []
        #integrate
        while s.successful() and s.t < tend:
                s.integrate(s.t + dt)
                solnt.append(s.t)
                solny.append(s.y/np.sum(s.y*np.conj(s.y))**.5)
        #store
        Cstor[J,M,:] = np.transpose(solny)[:,-1]
        Cstor[J,-M,:] = np.transpose(solny)[:,-1]

## Expectation value, incoherent, thermal average.
for J in range(Jmax):
    for M in range(-J,J+1):
        for jj in range(Jmax-2):
            w = 4*jj+6
            phi = np.angle(Cstor[J,M,jj])-np.angle(Cstor[J,M,jj+2])
            cos2 += Jweight[J]/(2*J+1)*(abs(Cstor[J,M,jj])**2*c2(jj,M) +
                    abs(Cstor[J,M,jj])*abs(Cstor[J,M,jj+2])*cp2(jj,M)*np.cos(w*tt+phi))

## End program
elapsed = time.time() - timer
print('\n Program took ' + str(round(elapsed)) + ' s to run. \n')
print('\n' + molecule+' at '+ str(I) + ' x 10$^{14}$ W cm$^{-2}$ at ' + str(TemperatureK) + ' K\n')

#Plot result, <cos**2\theta>
plt.figure()
plt.plot(tt*hbar/B*10**12,np.real(cos2),'k-')
plt.xlabel('Time [ps]')
plt.ylabel('<cos$^2\Theta$>')
plt.title(molecule+' at '+ str(I) + ' x 10$^{14}$ W cm$^{-2}$ at ' + str(TemperatureK) + ' K')
plt.grid()
plt.ylim(0,1)
plt.show()

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

# My libraries, if using ipython, restart if these are modified.
from boltzmann import Boltzmann
from pulse import pulse
from cosfun import *
from molecules import *
from const import *

# close old plots.
plt.close('all')

# start time
timer = time.time()

## tuneable parameters
molecule = 'N2'
pulse_FWHM = 200e-15 #FWHM duration of the intesity of the pulse in seconds
I = .1 #in 10**14 W cm**-2
TemperatureK = 70  #in Kelvin

##Calculated Parameters
#molecular parameters
B = B_dict[molecule]*1.98648e-23 #rotational constant in ground state in Joules
D = D_dict[molecule]*1.98648e-23 #centrifugal distorion in ground state in Joules
delta_alpha =  d_alpha_dict[molecule] #anisotropic polarizability
Jmax = Jmax_dict[molecule] #approximate max J
Temperature = k*TemperatureK/B
Jweight = Boltzmann(Temperature, 70, molecule) #Boltzmann distribution

#laser parameters
sigma = pulse_FWHM*B/hbar
E0 = 2.74*10**10*I**.5 # electric field amplitude
strength=0.5*4*np.pi*epsilon0*delta_alpha*E0**2/B #in rotational constan

## RHS of the Schrodinger Equation
def rhs(t, x, Mparm):
    dx = np.array(zeros(Jmax, dtype = 'complex'))
    Delta_omega = pulse(t, strength, sigma)
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
tend = 3*sigma; dt = .004*sigma
tt = np.linspace(0,8,1000)
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
            cos2 += Jweight[J]/(2*J+1)*(abs(Cstor[J,M,jj])**2*c2(jj,M) +abs(Cstor[J,M,jj])*abs(Cstor[J,M,jj+2])*cp2(jj,M)*np.cos(w*tt+phi))

## End program
elapsed = time.time() - timer
print('\n Program took ' + str(round(elapsed)) + ' s to run. \n')
print('\n' + molecule+' at '+ str(I) + ' x 10$^{14}$ W cm$^{-2}$ at ' + str(TemperatureK) + ' K\n')

#Plot result, <cos**2\theta>
plt.style.use('default')
plt.figure()
plt.plot(tt*hbar/B*10**12,np.real(cos2),'k-')
plt.xlabel('Time [ps]')
plt.ylabel('<cos$^2\Theta$>')
plt.title(molecule+' at '+ str(I) + ' x 10$^{14}$ W cm$^{-2}$ at ' + str(TemperatureK) + ' K')
plt.grid()
plt.ylim(0,1)
plt.show()



from pulse import pulse
ttt=np.linspace(-0.2*sigma,4*sigma,200)
ppp=np.zeros(np.size(ttt))
for ii in range(np.size(ttt)):
    ppp[ii]=pulse(ttt[ii],strength,sigma)
    
plt.plot(ttt,ppp)

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from SI2au import SI2au


import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import ode
import time
import math
import seaborn as sns
#import SI2au as SI2au

# My libraries, if using ipython, restart if these are modified.
from boltzmann import Boltzmann
from cosfun import *
from molecules import *
from const import *
from SI2au import SI2au

sns.set();


# todo显示范围与时间间隔
#注意强度E0未


def stov(tnc,ync,strength=98.25760788395505,TOAM=1,lambdaz=400e-9,ttm=0,Bm=3.9530952e-23):
    
    
    
    c=2.9979248e8;
    #c=SI2au(c,'v')# a.u.
    T=lambdaz/c;
    #epsilon0=SI2au(8.84e-9,'epsilon');#a.u.
    nu0=1/(lambdaz/c); #a.u
    omega0=nu0*2*np.pi;
    k0=2*np.pi/lambdaz;
    w_yi=150*lambdaz;#束腰半径
    
    
    
    

    z=0;
    #lambdaz = SI2au(lambda_a, 'a0')#a.u.波长
    y=ync/lambdaz
    
    ttra=(tnc*hbar/Bm-ttm)/T;

    w_y = w_yi;
    Ep=0;phi_st=0
    yi = c*ttra*T - z*lambdaz;#纵向坐标
    
    
    ## 时空平面内相位
    if (yi == 0) and (y== 0):
        phi_st = 0;
    elif (yi == 0) and (y > 0):
        phi_st = np.pi/2;
    elif (yi == 0) and (y < 0):
        phi_st = -np.pi/2 + 2*np.pi;
    elif (yi < 0) and (y>= 0):
        phi_st =np.arctan((y*lambdaz*w_yi)/(yi*w_y)) + np.pi;
    elif (yi < 0) and (y < 0):
        phi_st = np.arctan((y*lambdaz*w_yi)/(yi*w_y)) - np.pi + 2*np.pi;
    elif (yi > 0) and (y< 0):
        phi_st = np.arctan((y*lambdaz*w_yi)/(yi*w_y)) + 2*np.pi;
    else:
        phi_st = np.arctan((y*lambdaz*w_yi)/(yi*w_y));
        
    pp =strength*((yi/w_yi)**2 + (y*lambdaz/w_y)**2)**(np.abs(TOAM)/2)*np.exp((-1)*((yi/w_yi)**2 + ((y*lambdaz)/w_y)**2))*np.exp(1j*((-1)*TOAM*phi_st+omega0*ttra*T+ k0*z*lambdaz));

    
    p_amplitude = np.abs(pp);
    p_phi = np.angle(pp);
    p=p_amplitude*np.exp(1j*p_phi);
    phi_st=-1*TOAM*phi_st
    phi_st=phi_st%(2*np.pi);

    if phi_st > np.pi and phi_st<= 2*np.pi:
        phi_st= phi_st- 2*np.pi;
    elif phi_st< -np.pi and phi_st>= -2*np.pi:
        phi_st = phi_st+ 2*np.pi;

    p_phi = (p_phi- omega0 * ttra* T)%(2*np.pi);

    if p_phi > np.pi and p_phi<= 2*np.pi:
        p_phi = p_phi - 2*np.pi;
    elif p_phi< -np.pi and p_phi >= -2*np.pi:
        p_phi= p_phi+ 2*np.pi;

    
    return p,p_amplitude,p_phi,phi_st
    
    
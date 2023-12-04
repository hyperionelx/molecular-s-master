#import phi_st as phi_st
#import SI2au as SI2au

# def pulse(t, strength, sigma):
#     if (t >= 0  and  t <= 2*sigma):
#         p = strength*np.sin(np.pi*t/2/sigma)**2
#     else:
#         p = 0
#     return p
# 特定时间点不同位置处对应的方位角: 前提是认为 STOV 在 x 轴和 t 轴上同样宽（标准的圆环）

# def myfunc_phi(x, t):# ...之前转换的 myfunc_phi 函数的内容...


import numpy as np

#import library
#from SI2au import *
from phi_st import phi_st


def EA(x, t, c, lambda_au, omega, TOAM):
    phi=phi_st(x,t)
    # 方位角
    correct = 0.278324362050892
    envelope = np.sqrt((c**2*t**2/(4*lambda_au)**2)+x**2/(4*lambda_au)**2)*np.exp(-c**2*t**2/(4*lambda_au)**2-x**2/(4*lambda_au)**2)
    TV = correct*envelope
    if (t >= 0  and  t <= 2*sigma):
        Ex = correct*envelope*np.cos((-1)*TOAM*phi+omega*t)
        Ey = 0
        Ax = -TV/omega*np.sin((-1)*TOAM*phi+omega*t)
        Ay = 0

    return Ex, Ey, Ax, Ay



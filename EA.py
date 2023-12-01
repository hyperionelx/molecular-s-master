import numpy as np
import phi_st as phi_st
import SI2au as SI2au

# def pulse(t, strength, sigma):
#     if (t >= 0  and  t <= 2*sigma):
#         p = strength*np.sin(np.pi*t/2/sigma)**2
#     else:
#         p = 0
#     return p
# 特定时间点不同位置处对应的方位角: 前提是认为 STOV 在 x 轴和 t 轴上同样宽（标准的圆环）

# def myfunc_phi(x, t):# ...之前转换的 myfunc_phi 函数的内容...



def EA(x, t, c, lambda_au, omega, TOAM):
    lambda_m = 400e-9  # 波长，单位：米# 定义 lambda
    lambda_au = SI2au(lambda_m, 'a0')# 调用 SI2au 函数进行单位转换

    # 方位角
    correct = 0.278324362050892
    envelope = np.sqrt((c**2 * t**2 / (4 * lambda_au)**2) + x**2 / (4 * lambda_au)**2) * np.exp(-c**2 * t**2 / (4 * lambda_au)**2 - x**2 / (4 * lambda_au)**2)
    TV = correct * envelope
    Ex = TV * np.cos(-TOAM * phi_st + omega * t)
    Ey = 0
    Ax = -TV / omega * np.sin(-TOAM * phi_st + omega * t)
    Ay = 0

    return Ex, Ey, Ax, Ay



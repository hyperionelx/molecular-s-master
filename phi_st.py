import numpy as np

def phi_st(x, t):
    phi_st = np.zeros(len(x))

    for j in range(len(x)):
        # 时空平面内的螺旋相位
        if t == 0 and x[j] == 0:
            phi_st[j] = 0
        elif t == 0 and x[j] > 0:
            phi_st[j] = np.pi / 2
        elif t == 0 and x[j] < 0:
            phi_st[j] = -np.pi / 2 + 2 * np.pi
        elif t < 0 and x[j] >= 0:
            phi_st[j] = np.arctan(x[j] / t) + np.pi
        elif t < 0 and x[j] < 0:
            phi_st[j] = np.arctan(x[j] / t) - np.pi + 2 * np.pi
        elif t > 0 and x[j] < 0:
            phi_st[j] = np.arctan(x[j] / t) + 2 * np.pi
        else:
            phi_st[j] = np.arctan(x[j] / t)

    return phi_st


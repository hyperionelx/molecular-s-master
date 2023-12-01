import numpy as np

def pulse(t, strength, sigma):
    if (t >= 0  and  t <= 2*sigma):
        p = strength*np.sin(np.pi*t/2/sigma)**2
    else:
        p = 0
    return p

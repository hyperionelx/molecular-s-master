##Wigner 3 J symbol using Racah Formula
#Craig Benko 2014.04.30
from pylab import *
import math as m
import seaborn as sns
import numpy as np


def Wigner3J(j1, j2, j3, m1, m2, m3):
    if 2*j1 != floor(2*j1) or 2*j2 != floor(2*j2) or 2*j3 != floor(2*j3) or 2*m1 != floor(2*m1) or 2*m2 != floor(2*m2) or 2*m3 != floor(2*m3):
        # print('All arguments must be integers')
        wig = 0
    elif j1 - m1 != floor(j1 -m1):
        # print('2*j1 must have same parity as 2*m1')
        wig = 0
    elif j2 - m2 != floor(j2 -m2):
        # print('2*j2 must have same parity as 2*m2')
        wig = 0
    elif j3 - m3 != floor(j3 -m3):
        # print('2*j3 must have same parity as 2*m3')
        wig = 0
    elif j3 > j1 + j2 or j3 < abs(j1 - j2):
        # print('j3 out of bounds')
        wig = 0
    elif abs(m1) > j1:
        # print('m1 is out of bounds')
        wig = 0
    elif abs(m2) > j2:
        # print('m2 is out of bounds')
        wig = 0
    elif abs(m3) > j3:
        # print('m3 is out of bounds')
        wig = 0
    else:
        t1 = j2 - m1 - j3
        t2 = j1 + m2 - j3
        t3 = j1 + j2 - j3
        t4 = j1 - m1
        t5 = j2 + m2

        tmin = max(0, max(t1, t2))
        tmax = min(t3, min(t4, t5))

        wig = 0

        for t in range(tmin, tmax+1):
            wig += (-1)**t / (m.factorial(t)*m.factorial(t-t1)*\
                m.factorial(t-t2)*m.factorial(t3-t)*m.factorial(t4-t)*\
                m.factorial(t5-t) )

        wig *= (-1)**(j1-j2-m3)*sqrt(m.factorial(j1+j2-j3)*\
            m.factorial(j1-j2+j3)*m.factorial(-j1+j2+j3)/\
            m.factorial(j1+j2+j3+1)*m.factorial(j1+m1)*\
            m.factorial(j1-m1)*m.factorial(j2+m2)*m.factorial(j2-m2)*\
            m.factorial(j3+m3)*m.factorial(j3-m3))

    return wig

from numpy import zeros, shape
## Normal
def c2(J,M):
    if abs(M)>J:
        return 0
    else:
        return 1/3 + 2/3*((J*(J+1)-3*M**2))/((2*J+3)*(2*J-1))

def cp2(J,M):
    if abs(M)>J:
        return 0
    else:
        return ((2.0*J+1)*(2*J+5)*(J+1-M))**.5*((2+J-M)*(J+1+M)*(J+2+M))**.5/((2*J+1)*(2*J+3)*(2*J+5))

def cm2(J,M):
    if abs(M)>J:
        return 0
    elif M==J:
        return 0
    else:
        return ((2.0*J-3)*(2*J+1)*(J-1-M))**.5*((J-M)*(J-1+M)*(J+M))**.5/((2*J-3)*(2*J-1)*(2*J+1))

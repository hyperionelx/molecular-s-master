import numpy as np

def Boltzmann(T, Jm, mol):
    # Note: the temperature needs to be in units of the rotational constant.
    weven_dict = {'O2': 0, 'N2': 2, 'N2O': 1, 'CO2': 1, 'D2': 2, 'OCS': 1}
    wodd_dict =  {'O2': 1, 'N2': 1, 'N2O': 1, 'CO2': 0, 'D2': 1, 'OCS': 1}

    weven = weven_dict[mol]
    wodd = wodd_dict[mol]
    pop = np.zeros(Jm)
    for J in range(Jm):
        if np.mod(J,2) == 0:
            pop[J] = weven*(2*J+1)*np.exp(-J*(J+1)/T)
        else:
            pop[J] = wodd*(2*J+1)*np.exp(-J*(J+1)/T)

    return pop/np.sum(pop)

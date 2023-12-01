# Field-free molecular alignment with Python

These files contain the basic constituents to simulate field-free molecular alignment from a single femtosecond pulse.  The simluation integrates the Schrodinger equation and performs an incoherent, thermal average to generate a <cos**2 \theta>(t) expectation value.

The code is basically pure Python and makes use of Scipy and Numpy.  The ODE solver is from Scipy.  To run the simulation, simply call 'align.py' from the command line.  Various parameters can be changed in the 'align.py' file.  For example, molecule, gas temperature, and laser parameters.

The code is currently pretty slow, it takes about 10 seconds to run the program for N2 (Jmax ~24).  For larger molecules like N2O, it takes a few minutes (Jmax ~60).

Future implementations will use Cython to speed up the code, but it is currently tolerable for N2.  I will also create a parameter file to keep all the commonly changed parameters currently located in 'align.py'.  I will also include calculations of more expecation values and the rotational state density.  


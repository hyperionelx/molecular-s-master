import numpy as np
import matplotlib.pyplot as plt

# Define the parameters
delta_alpha = 1.0  # The polarizability difference, arbitrary units
theta = np.pi/2    # Alignment angle, for maximal alignment use pi/2
Phi = np.pi/4      # Phase, arbitrary choice
tau = 1.0          # Time delay for the probe field, arbitrary units
omega_0 = 1.0      # Angular frequency, arbitrary units

# Time array
t = np.linspace(0, 10, 1000)  # Arbitrary time units

# Define the envelope function epsilon_0 as a function of (t - tau)
# For simplicity, assuming a Gaussian form for the envelope function
epsilon_0 = np.exp(-(t - tau)**2)

# Define the probe field E_pr(t)
E_pr = epsilon_0 * np.exp(-1j * omega_0 * (t - tau))

# Calculate the signal field E_sig(t) based on the given formula
E_sig = (E_pr / 4) * delta_alpha * np.sin(theta)**2 * np.exp(1j * 2 * Phi)

# Plotting the real and imaginary parts of the signal field
plt.figure(figsize=(12, 6))

# Real part of the signal
plt.subplot(1, 2, 1)
plt.plot(t, E_sig.real, label='Real part')
plt.title('Real part of the Signal Field')
plt.xlabel('Time (arbitrary units)')
plt.ylabel('Electric Field (arbitrary units)')
plt.legend()

# Imaginary part of the signal
plt.subplot(1, 2, 2)
plt.plot(t, E_sig.imag, label='Imaginary part')
plt.title('Imaginary part of the Signal Field')
plt.xlabel('Time (arbitrary units)')
plt.ylabel('Electric Field (arbitrary units)')
plt.legend()

plt.tight_layout()
plt.show()







# Updated parameters for the new conditions
t_fixed = 10.68e-12  # Fixed time set to 10.68 ps
tau_values = np.linspace(0, 100, 1000)  # Tau range from 0 to 100 (arbitrary units)

# Assuming the units for tau are the same as for t, if not, an appropriate conversion factor should be applied
# Reinitialize the array to store the values of E_sig at t=10.68 ps for different tau
E_sig_at_t_fixed = np.empty(len(tau_values), dtype=complex)

# Compute E_sig(t_fixed) for each value of tau
for i, tau in enumerate(tau_values):
    epsilon_0 = np.exp(-(t_fixed - tau)**2)  # Assuming Gaussian for epsilon_0
    E_pr_at_t_fixed = epsilon_0 * np.exp(-1j * omega_0 * (t_fixed - tau))
    E_sig_at_t_fixed[i] = (E_pr_at_t_fixed / 4) * delta_alpha * np.sin(theta)**2 * np.exp(1j * 2 * Phi)

# Plotting the real and imaginary parts of E_sig(t_fixed) as a function of tau
plt.figure(figsize=(12, 6))

# Real part of E_sig
plt.subplot(1, 2, 1)
plt.plot(tau_values, E_sig_at_t_fixed.real, label='Real part')
plt.title('Real part of E_sig(t_fixed) vs. Tau')
plt.xlabel('Tau (arbitrary units)')
plt.ylabel('E_sig(t_fixed) Real (arbitrary units)')
plt.legend()

# Imaginary part of E_sig
plt.subplot(1, 2, 2)
plt.plot(tau_values, E_sig_at_t_fixed.imag, label='Imaginary part')
plt.title('Imaginary part of E_sig(t_fixed) vs. Tau')
plt.xlabel('Tau (arbitrary units)')
plt.ylabel('E_sig(t_fixed) Imaginary (arbitrary units)')
plt.legend()

plt.tight_layout()
plt.show()

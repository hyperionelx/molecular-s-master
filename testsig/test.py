import numpy as np
import matplotlib.pyplot as plt


# Define the parameters
delta_alpha = 1.0  # The polarizability difference, arbitrary units
omega_0 = 1.0      # Angular frequency, arbitrary units
theta = np.pi/2    # Alignment angle, for maximal alignment use pi/2
Phi = np.pi/4      # Phase, arbitrary choice

# Updated parameters for the new conditions
t_fixed = 10.68e-12  # Fixed time set to 10.68 ps
tau_fixed = 50e-12
tau_values = np.linspace(0, 100, 1000) * 1e-12  # Tau range from 0 to 100 (arbitrary units)
# Define the time array for the range of interest
t_values = np.linspace(0, 100, 1000) * 1e-12  # Convert range to picoseconds

# Assuming the units for tau are the same as for t, if not, an appropriate conversion factor should be applied
# Reinitialize the array to store the values of E_sig at t=10.68 ps for different tau
E_sig_at_t_fixed = np.empty(len(tau_values), dtype=complex)
E_pr_at_t_fixed = np.empty(len(t_values), dtype=complex)


# Parameters
tau = 10.68e-12  # Fixed tau value in seconds
omega_0 = 2 * np.pi * 1e12  # Angular frequency in radians per second, assuming 1 THz

# Define the time domain
t = np.linspace(0, 5e-11, 1000)  # Time from 0 to 50 ps

# Define the Gaussian function for epsilon_0
def epsilon_0(t, tau, width=1e-12):
    return np.exp(-((t - tau) ** 2) / (2 * width ** 2))

# Calculate E_pr(t)
E_pr = epsilon_0(t, tau) * np.exp(-1j * omega_0 * (t - tau))
E_pr_at_t_fixed = epsilon_0(t_fixed, tau) * np.exp(-1j * omega_0 * (t_fixed - tau))




# Compute E_sig(t_fixed) for each value of tau
for i, tau in enumerate(tau_values):
    epsilon_0 = np.exp(-(t_fixed - tau)**2)  # Assuming Gaussian for epsilon_0
    E_pr_at_t_fixed = epsilon_0 * np.exp(-1j * omega_0 * (t_fixed - tau))
    E_sig_at_t_fixed[i] = (E_pr_at_t_fixed / 4) * delta_alpha * np.mean(np.sin(theta)**2 * np.exp(1j * 2 * Phi))

# Plotting the real and imaginary parts of E_sig(t_fixed) as a function of tau
plt.figure(figsize=(12, 6))

# # Real part of E_sig
# plt.subplot(1, 2, 1)
# plt.plot(tau_values, E_sig_at_t_fixed.real, label='Real part')
# plt.title('Real part of E_sig(t_fixed) vs. Tau')
# plt.xlabel('Tau (arbitrary units)')
# plt.ylabel('E_sig(t_fixed) Real (arbitrary units)')
# plt.legend()

# # Imaginary part of E_sig
# plt.subplot(1, 2, 2)
# plt.plot(tau_values, E_sig_at_t_fixed.imag, label='Imaginary part')
# plt.title('Imaginary part of E_sig(t_fixed) vs. Tau')
# plt.xlabel('Tau (arbitrary units)')
# plt.ylabel('E_sig(t_fixed) Imaginary (arbitrary units)')
# plt.legend()

# Total of E_sig
plt.plot(tau_values, E_sig_at_t_fixed, label='Imaginary part')
plt.title('Imaginary part of E_sig(t_fixed) vs. Tau')
plt.xlabel('Tau (arbitrary units)')
plt.ylabel('E_sig(t_fixed) Imaginary (arbitrary units)')
plt.legend()

plt.tight_layout()
plt.show()

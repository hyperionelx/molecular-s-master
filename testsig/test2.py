import numpy as np
import matplotlib.pyplot as plt


# Define the parameters
delta_alpha = 1.0  # The polarizability difference, arbitrary units
omega_0 = 1.0      # Angular frequency, arbitrary units
theta = np.pi/2    # Alignment angle, for maximal alignment use pi/2
Phi = np.pi/4      # Phase, arbitrary choice


# Define a range of tau values
tau_values = np.linspace(0, 100, 1)
# Define the time array for the range of interest
t_values = np.linspace(0, 100, 1000) * 1e-12  # Convert range to picoseconds
# t_values = 40e-12


# Define a new envelope function epsilon_0 as a Gaussian centered at tau
def epsilon_0(t, tau):
    # Assuming a width for the Gaussian pulse
    width = 1e-12  # Width of the Gaussian in picoseconds
    return np.exp(-((t - tau) ** 2) / (2 * width ** 2))

# Now we'll compute E_sig for a fixed phi and a range of tau
phi = np.pi / 4  # Phase, arbitrary choice

# Compute E_sig for each combination of t and tau
E_sig = np.empty((len(t_values), len(tau_values)), dtype=complex)

for j, t in enumerate(t_values):
    for i, tau in enumerate(tau_values):
        E_sig[j, i] = delta_alpha * epsilon_0(t, tau) * np.exp((-1.0j) * (omega_0 * (t - tau) - phi))

# For visualization, we'll plot the absolute value of E_sig at t=10.68 ps
t_fixed_index = np.argmin(np.abs(t_values - 10.68e-12))  # Find the index of the closest value to 10.68 ps

# Plotting the magnitude of E_sig(t_fixed) as a function of tau
plt.figure(figsize=(10, 5))
plt.plot(tau_values, np.abs(E_sig[t_fixed_index, :]), label='|E_sig(t_fixed)|')
plt.title('Magnitude of E_sig(t_fixed) vs. Tau')
plt.xlabel('Tau (arbitrary units)')
plt.ylabel('|E_sig(t_fixed)| (arbitrary units)')
plt.legend()
plt.show()

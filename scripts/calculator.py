import numpy as np

k_B = 1.380649e-23
c = 299792458
hbar = 1.054571817e-34
T_CMB = 2.735
zeta = 1.202056903159594285399738161511449990764986292
Mpc = 3.08567758e22
G = 6.6743e-11
Omega_b0 = 0.05
h = 0.7
km = 1000
mH = 1.6735575e-27

n_gam0 = 2/np.pi**2 * zeta * (T_CMB*k_B/(hbar*c))**3
print(f"{n_gam0:.3e}")

n_b0 = 3 * Omega_b0 * (100*h*km/Mpc)**2 / (8 * np.pi * G * mH)
print(f"{n_b0:.3e}")

print(f"{n_b0/n_gam0:.3e}")
print(f"{n_b0/n_gam0/Omega_b0/h**2:.3e}")
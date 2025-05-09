import numpy as np
import healpy as hp
import matplotlib.pyplot as plt
plt.rcParams.update({'text.usetex': True, 'font.size': 18, 'font.family': 'serif', 'font.serif': 'Computer Modern Sans Serif', 'font.weight': 100, 'mathtext.fontset': 'cm', 'xtick.labelsize': 16, 'ytick.labelsize': 16})

# Load data
data = np.loadtxt("results/C_ells.txt")
ells = data[:, 0]
D_ell = data[:, 1]

# Create full C_ell array including ell = 0, 1 (set to zero)
ell_max = int(ells[-1])
C_ell = np.zeros(ell_max + 1)
C_ell[2:] = 2*np.pi / (ells*(ells + 1)) * D_ell

# Generate a map from C_ell
nside = 4096*2 
CMB_map = hp.synfast(C_ell, nside = nside, lmax = ell_max, new = True)

# Plot the map
rms = np.std(CMB_map)
hp.mollview(CMB_map, title = "Simulated CMB Map", min = -3*rms, max = 3*rms, unit = r"$\mu$K", cmap = "nipy_spectral")
plt.show()

hp.mollview(CMB_map, title = "Simulated CMB Map")
plt.show()

# Save the map as fits file
hp.write_map("results/CMB_map.fits", CMB_map, overwrite = True)


# # Cl is a 3 x (lmax+1) array: rows are TT, EE, BB, TE
# Cl = [C_TT, C_EE, C_BB, C_TE]  # Add zeros if you don't use all

# nside = 512
# lmax = len(C_TT) - 1

# # Generate T, Q, U maps
# T_map, Q_map, U_map = hp.synfast(Cl, nside=nside, lmax=lmax, new=True, pol=True)

# # Assuming you already have Q_map and U_map
# maps = [np.zeros_like(Q_map), Q_map, U_map]  # T, Q, U
# alm = hp.map2alm(maps, pol=True)

# # Convert to E and B alm
# a_T, a_E, a_B = alm  # These are the spherical harmonic coefficients

# # Create E-mode and B-mode maps (set T=0 for visualization)
# zero_T = np.zeros_like(a_E)
# map_E = hp.alm2map([zero_T, a_E, np.zeros_like(a_E)], nside=hp.get_nside(Q_map), pol=True)[1]
# map_B = hp.alm2map([zero_T, np.zeros_like(a_B), a_B], nside=hp.get_nside(Q_map), pol=True)[2]

# hp.mollview(map_E, title='E-mode Polarization', unit='μK', norm='hist', cmap='coolwarm')
# plt.show()

# hp.mollview(map_B, title='B-mode Polarization', unit='μK', norm='hist', cmap='coolwarm')
# plt.show()
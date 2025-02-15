import numpy as np
import matplotlib.pyplot as plt

"""
Calculate analytical expressions
"""
c = 299792458               # [m/s]
Mpc = 3.08567758e22         # [m]
Myr = 3.1536e13             # [s]
Gyr = 3.1536e16             # [s]

# Hubble parameter today
h = 0.67
H0 = 100*h                  # [km/sec/Mpc]

# Density parameters today
Omega_m0 = 0.05 + 0.267
Omega_r0 = 5.50896e-05 + 3.81093e-05
Omega_Lambda0 = 0.682907

# Radiation-matter equality
x_rm = np.log(Omega_r0/Omega_m0)
z_rm = Omega_m0/Omega_r0 - 1
t_rm = Omega_r0**(3/2) / (2*H0*Omega_m0**2) * (Mpc/1000/Myr)
print("Radiation-matter equality:")
print(f"x: {x_rm:.3f}")
print(f"z: {z_rm:.3f}")
print(f"t: {t_rm:.3f} Myr\n")

# Matter-dark energy equality
x_mLambda = np.log(Omega_m0/Omega_Lambda0) / 3
z_mLambda = (Omega_Lambda0/Omega_m0)**(1/3) - 1
t_mLambda = 1/(3*H0) * (2/np.sqrt(Omega_Lambda0) - Omega_r0**(3/2) / (2*Omega_m0**2)) * (Mpc/1000/Gyr)
print("Matter-dark energy equality:")
print(f"x: {x_mLambda:.3f}")
print(f"z: {z_mLambda:.3f}")
print(f"t: {t_mLambda:.3f} Gyr\n")

# Onset of acceleration
x_acc = np.log(Omega_m0/(2*Omega_Lambda0)) / 3
z_acc = (2*Omega_Lambda0/Omega_m0)**(1/3) - 1
t_acc = 1/(3*H0) * (np.sqrt(2/Omega_Lambda0) - Omega_r0**(3/2) / (2*Omega_m0**2)) * (Mpc/1000/Gyr)
print("Onset of acceleration:")
print(f"x: {x_acc:.3f}")
print(f"z: {z_acc:.3f}")
print(f"t: {t_acc:.3f} Gyr\n")


"""
Read data from file
"""
n_pts = 100_000
data = {}
keys = ["x", "eta", "Hp", "dHpdx", "ddHpddx", "Omega_b", "Omega_CDM", "Omega_Lambda", "Omega_gamma", "Omega_nu", "Omega_k", "t", "detadx", "chi", "d_L", "d_A", "T_CMB"]
for j, key in enumerate(keys):
    data[key] = np.zeros(n_pts)

with open("results/cosmology.txt", "r") as infile:
    for i, line in enumerate(infile):
        values = line.split()
        for j, key in enumerate(keys):
            data[key][i] = float(values[j])


"""
Plot of eta'Hp/c to test for numerical stability
"""
plt.plot(data["x"], data["detadx"]*data["Hp"]/c, "k")
plt.plot([data["x"][0], data["x"][-1]], [1, 1], color = "#de82b4", linestyle = "--")
plt.show()


"""
Plot Hp'/Hp and Hp''/Hp vs. x together with analytical expectations
"""
#TODO fix subplots
plt.plot(data["x"], data["dHpdx"]/data["Hp"], "k")
plt.plot([data["x"][0], x_rm, x_rm, x_mLambda, x_mLambda, data["x"][-1]], [-1, -1, -1/2, -1/2, 1, 1], color = "#de82b4", linestyle = "--")
# plt.plot([data["x"][0], x_rm], [-1, -1], color = "#de82b4")
# plt.plot([x_rm, x_mLambda], [-1/2, -1/2], color = "#de82b4")
# plt.plot([x_mLambda, data["x"][-1]], [1, 1], color = "#de82b4")
# plt.plot(data["x"], data["dHpdx"]/data["Hp"], "k")
plt.show()

#TODO factor 2 wrong in expression for double derivative?
plt.plot(data["x"], data["ddHpddx"]/data["Hp"], "k")
plt.plot([data["x"][0], x_rm, x_rm, x_mLambda, x_mLambda, data["x"][-1]], [1, 1, 1/4, 1/4, 1, 1], color = "#de82b4", linestyle = "--")
# plt.plot([data["x"][0], x_rm], [1, 1], color = "#de82b4")
# plt.plot([x_rm, x_mLambda], [1/4, 1/4], color = "#de82b4")
# plt.plot([x_mLambda, data["x"][-1]], [1, 1], color = "#de82b4")
# plt.plot(data["x"], data["dHpdx"]/data["Hp"], "k")
plt.show()


"""
Plot a vs. t together with analytical expectations
"""
#TODO something wrong here. how visualize?
#TODO argue that it is inaccurate in the matter dominated era since radiation is not negligible in the beginning and dark energy not in the end?
t_r = np.linspace(data["t"][0], t_rm * Myr, 100) / Gyr
t_m = np.linspace(t_rm / Gyr, t_mLambda, 100)
t_Lambda = np.linspace(t_mLambda, data["t"][-1] / Gyr, 100)

x_r = np.exp(np.linspace(data["x"][0], x_rm, 100))
x_m = np.exp(np.linspace(x_rm, x_mLambda, 100))
x_Lambda = np.exp(np.linspace(x_mLambda, data["x"][-1], 100))

plt.plot(data["t"] / Gyr, np.exp(data["x"]), "k")
plt.plot(t_r, x_r, color = "#de82b4", linestyle = "--")
plt.plot(t_m, x_m, color = "#de82b4", linestyle = "--")
plt.plot(t_Lambda, x_Lambda, color = "#de82b4", linestyle = "--")
plt.xscale("log")
plt.yscale("log")
plt.show()
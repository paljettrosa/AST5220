import numpy as np
import matplotlib.pyplot as plt

Gpc = 3.08567758e25

x = []
d_L = []
with open("results/cosmology.txt", "r") as infile:
    for line in infile:
        values = line.split()
        x.append(float(values[0]))
        d_L.append(float(values[-3])/Gpc)
z = 1/np.exp(np.array(x)) - 1
d_L = np.array(d_L)

z_measure = []
d_L_measure = []
d_L_err = []
with open("data/supernovadata.txt", "r") as infile:
    infile.readline()
    for line in infile:
        values = line.split()
        z_measure.append(float(values[0]))
        d_L_measure.append(float(values[1]))
        d_L_err.append(float(values[2]))
z_measure = np.array(z_measure)
d_L_measure = np.array(d_L_measure)
d_L_err = np.array(d_L_err)

plt.plot(z, d_L)
plt.errorbar(z_measure, d_L_measure, d_L_err)
plt.xlim(z_measure[0], z_measure[-1])
plt.ylim(min(d_L_measure-d_L_err), max(d_L_measure+d_L_err))
plt.show()


"""
Minimum chi^2 found 29.2811 0.70189 0.25932 0.0673887
"""
chi2 = []
h = []
Omega_m = []
Omega_k = []
with open("results/results_supernovafitting.txt", "r") as infile:
    infile.readlines(201)
    for line in infile:
        values = line.split()
        chi2.append(float(values[0]))
        h.append(float(values[1]))
        Omega_m.append(float(values[2]))
        Omega_k.append(float(values[3]))
chi2 = np.array(chi2)
H0 = 100*np.array(h)
Omega_m = np.array(Omega_m)
Omega_k = np.array(Omega_k)
Omega_Lambda = 1 - Omega_m - Omega_k

#TODO we have four parameters (H, Omega_m, Omega_k and Omega_Lambda), and one constraint (Omega_m + Omega_k + Omega_Lambda = 1), hence k=3 degrees of freedom
sigma2 = 8.02
sigma1 = 3.53
plt.scatter(Omega_m[chi2 < (np.min(chi2) + sigma2)], Omega_Lambda[chi2 < (np.min(chi2) + sigma2)])
plt.scatter(Omega_m[chi2 < (np.min(chi2) + sigma1)], Omega_Lambda[chi2 < (np.min(chi2) + sigma1)])
plt.scatter(0.25932, 1 - 0.25932 - 0.0673887)
plt.plot(Omega_m, 1 - Omega_m, "k--")
plt.show()

#TODO make subplots
std_H0 = np.std(H0)
print(std_H0)
plt.hist(H0[chi2 < np.min(chi2) + sigma2], bins = 50, color = "#50b8e7", edgecolor = "black")
plt.hist(H0[chi2 < np.min(chi2) + sigma1], bins = 30, color = "#b9e2f5", edgecolor = "black")
plt.axvline(H0[np.argmin(chi2)], color = "k", linestyle = "--")
plt.show()


std_Omega_m = np.std(Omega_m)
print(std_Omega_m)
plt.hist(Omega_m[chi2 < np.min(chi2) + sigma2], bins = 50, color = "#d84528", edgecolor = "black")
plt.hist(Omega_m[chi2 < np.min(chi2) + sigma1], bins = 30, color = "#eb8e27", edgecolor = "black")
plt.axvline(Omega_m[np.argmin(chi2)], color = "k", linestyle = "--")
plt.show()


std_Omega_k = np.std(Omega_k)
print(std_Omega_k)
plt.hist(Omega_k[chi2 < np.min(chi2) + sigma2], bins = 50, color = "#c94475", edgecolor = "black")
plt.hist(Omega_k[chi2 < np.min(chi2) + sigma1], bins = 30, color = "#e691b2", edgecolor = "black")
plt.axvline(Omega_k[np.argmin(chi2)], color = "k", linestyle = "--")
plt.show()


std_Omega_Lambda = np.std(Omega_Lambda)
print(std_Omega_Lambda)
plt.hist(Omega_Lambda[chi2 < np.min(chi2) + sigma2], bins = 50, color = "#9849b3", edgecolor = "black")
plt.hist(Omega_Lambda[chi2 < np.min(chi2) + sigma1], bins = 30, color = "#cd96e0", edgecolor = "black")
plt.axvline(1 - Omega_m[np.argmin(chi2)] - Omega_k[np.argmin(chi2)], color = "k", linestyle = "--")
plt.show()
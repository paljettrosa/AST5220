import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe
plt.rcParams.update({'text.usetex': True, 'font.size': 18, 'font.family': 'serif', 'font.serif': 'Computer Modern Sans Serif', 'font.weight': 100, 'mathtext.fontset': 'cm', 'xtick.labelsize': 16, 'ytick.labelsize': 16})

Gpc = 3.08567758e25
Omega_r0 = 5.50896e-05 + 3.81093e-05

"""
Read data from files
"""
x = []
d_L = []
with open("results/cosmology.txt", "r") as infile:
    for line in infile:
        values = line.split()
        x.append(float(values[0]))
        d_L.append(float(values[-3])/Gpc)
z = 1/np.exp(np.array(x)) - 1
d_L = np.array(d_L)

d_L_bestfit = []
with open("results/cosmology_bestfit.txt", "r") as infile:
    for line in infile:
        values = line.split()
        d_L_bestfit.append(float(values[-2])/Gpc)
d_L_bestfit = np.array(d_L_bestfit)

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
Omega_Lambda = 1 - Omega_m - Omega_k - Omega_r0


"""
Plot of luminosity distance versus redshift
"""
plt.figure(figsize = (8, 6))
plt.subplots_adjust(left = 0.09, right = 0.96, top = 0.93)
plt.plot(z, d_L_bestfit/z, "k", label = r"Minimum $\chi^2$")
plt.plot(z, d_L/z, "slategrey", label = "Planck")
plt.errorbar(z_measure, d_L_measure/z_measure, d_L_err/z_measure, fmt = ' ', elinewidth = 1.5, capsize = 3, ecolor = "#ff77bc", label = "Measurements")
plt.legend()
plt.xscale("log")
plt.xlim(10**(1.18*np.log10(z_measure[-1])), 10**(1.01*np.log10(z_measure[0])))
plt.ylim(0.99*np.min(d_L_measure/z_measure-d_L_err/z_measure), 1.005*np.max(d_L_measure/z_measure+d_L_err/z_measure))
plt.xticks([1, 0.7, 0.3, 0.1, 0.07, 0.03], labels = [1, 0.7, 0.3, 0.1, 0.07, 0.03])
plt.xlabel(r"Redshift $z$")
plt.ylabel(r"Luminosity distance $d_L/z$ [Gpc]")
plt.savefig("figs/luminosity_distance.pdf")
plt.show()


"""
Plot of the MCMC fits in the (Omega_Lambda, Omega_m)-plane
"""
sigma2 = 8.02
sigma1 = 3.53
plt.show()
plt.figure(figsize = (8, 6))
plt.subplots_adjust(left = 0.09, right = 0.96, top = 0.93)
plt.scatter(Omega_m[chi2 < (np.min(chi2) + sigma2)], Omega_Lambda[chi2 < (np.min(chi2) + sigma2)], s = 15, color = "#de82b4", label = r"$\chi^2<2\sigma$ fits")
plt.scatter(Omega_m[chi2 < (np.min(chi2) + sigma1)], Omega_Lambda[chi2 < (np.min(chi2) + sigma1)], s = 15, color = "#8d3063", label = r"$\chi^2<1\sigma$ fits")
plt.scatter(Omega_m[np.argmin(chi2)], Omega_Lambda[np.argmin(chi2)], color = "black", label = r"Minimum $\chi^2$")
plt.scatter(0.317, 0.682, color = "slategrey", label = "Planck", zorder = 2)
plt.plot(np.linspace(0, 1, 100), 1 - np.linspace(0, 1, 100), "k--", label = "Flat universe", zorder = 1)
plt.legend(framealpha = 1)
plt.xlim(0, 1)
plt.xlabel(r"Matter density $\Omega_{m0}$")
plt.ylabel(r"Dark energy density $\Omega_{\Lambda0}$")
plt.savefig("figs/MCMC_fits.pdf")
plt.show()



"""
Plot of posterior PDFs for H0, Omega_m0, Omega_k0 and Omega_Lambda0
"""
variables = [[H0, Omega_m], [Omega_k, Omega_Lambda]]
Planck = [[67.4, 0.315], [0, 0.685]]
names = [["H0", "Omega_m0"], ["Omega_k0", "Omega_Lambda0"]]
labels = [[r"$H_0$ [km/s/Mpc]", r"$\Omega_{m0}$"], [r"$\Omega_{k0}$", r"$\Omega_{\Lambda0}$"]]
colors = [["#b9e2f5", "#f5b267"], ["#f2a2c1", "#ddabed"]]

fig, axs = plt.subplots(2, 2, figsize = (16, 12), constrained_layout = True)
for i in range(2):
    for j in range(2):
        mean = np.mean(variables[i][j][chi2 < np.min(chi2) + sigma1])
        std = np.std(variables[i][j][chi2 < np.min(chi2) + sigma1])
        print(f"  {names[i][j]}:")
        print(f"Mean:     {mean:6.2f} km/s/Mpc" if (i == 0 and j == 0) else f"Mean:     {mean:6.4f}")
        print(f"Std:      {std:6.2f} km/s/Mpc" if (i == 0 and j == 0) else f"Std:      {std:6.4f}")
        print(f"Best fit: {variables[i][j][np.argmin(chi2)]:6.2f} km/s/Mpc" if (i == 0 and j == 0) else f"Best fit: {variables[i][j][np.argmin(chi2)]:6.4f}")
        print(f"Planck:   {Planck[i][j]:6.2f} km/s/Mpc" if (i == 0 and j == 0) else f"Planck:   {Planck[i][j]:6.4f}")
        print("\n")

        x = np.linspace(mean - 3*std, mean + 3*std, 1000)
        pdf = 1/(np.sqrt(2*np.pi)*std) * np.exp(-(x-mean)**2/(2*std**2))

        axs[i][j].hist(variables[i][j][chi2 < np.min(chi2) + sigma1], density = True, bins = 30, color = colors[i][j], edgecolor = "black")
        axs[i][j].plot(x, pdf, color = "black", label = "Posterior fit")
        axs[i][j].axvline(variables[i][j][np.argmin(chi2)], color = "crimson", label = r"Minimum $\chi^2$", linewidth = 2, path_effects=[pe.Stroke(linewidth = 4, foreground = "black"), pe.Normal()])
        axs[i][j].axvline(mean, color = "forestgreen", label = r"Mean", linewidth = 2, path_effects=[pe.Stroke(linewidth = 4, foreground = "black"), pe.Normal()])
        axs[i][j].axvline(Planck[i][j], color = "dodgerblue", label = r"Planck best-fit", linewidth = 2, path_effects=[pe.Stroke(linewidth = 4, foreground = "black"), pe.Normal()])
        axs[i][j].set_xlabel(labels[i][j])

axs[0][0].legend(framealpha = 1)
fig.supylabel("Probability distribution functions")
fig.savefig("figs/distributions.pdf")
plt.show()
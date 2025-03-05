import numpy as np
import matplotlib.pyplot as plt
plt.rcParams.update({'text.usetex': True, 'font.size': 18, 'font.family': 'serif', 'font.serif': 'Computer Modern Sans Serif', 'font.weight': 100, 'mathtext.fontset': 'cm', 'xtick.labelsize': 16, 'ytick.labelsize': 16})

"""
Calculate analytical expressions
"""
c = 299792458               # [m/s]
Mpc = 3.08567758e22         # [m]
Gpc = 3.08567758e25         # [m]
yr = 3.1536e7               # [s]
kyr = 3.1536e10             # [s]
Myr = 3.1536e13             # [s]
Gyr = 3.1536e16             # [s]

# Hubble parameter today
h = 0.67
H0 = 100*h                  # [km/s/Mpc]

# Density parameters today
Omega_m0 = 0.05 + 0.267
Omega_r0 = 5.50896e-05 + 3.81093e-05
Omega_Lambda0 = 0.682907

"""
Read data from file
"""
x_min = -12
x_max = 0
n_pts = (x_max - x_min)*100 + 1
data = {}
keys = ["x", "Xe", "ne", "tau", "dtaudx", "ddtauddx", "g_tilde", "dgdx_tilde", "ddgddx_tilde"]
for j, key in enumerate(keys):
    data[key] = np.zeros(n_pts)

with open("results/recombination.txt", "r") as infile:
    for i, line in enumerate(infile):
        values = line.split()
        for j, key in enumerate(keys):
            data[key][i] = float(values[j])

plt.plot(data["x"], data["tau"])
plt.plot(data["x"], -data["dtaudx"])
plt.plot(data["x"], data["ddtauddx"])
plt.xlim(x_min, x_max)
plt.yscale("log")
plt.show()

plt.plot(data["x"], data["Xe"])
plt.yscale("log")
plt.show()

#TODO remember to test that the integral of g_tilde is 1!
plt.plot(data["x"], data["g_tilde"])
plt.show()

plt.plot(data["x"], data["dgdx_tilde"])
plt.show()

plt.plot(data["x"], data["ddgddx_tilde"])
plt.show()




# """
# Plot of eta'Hp/c to test for numerical stability
# """
# plt.figure(figsize = (8, 6))
# plt.subplots_adjust(left = 0.16, right = 0.95, top = 0.93)
# plt.scatter(data["x"], data["detadx"]*data["Hp"]/c, 15, color = "black")
# plt.plot([data["x"][0], data["x"][-1]], [1, 1], color = "#ff77bc")
# plt.xlim([x_min, x_max])
# plt.xlabel(r"$x=\log(a)$")
# plt.ylabel(r"$\frac{\eta'\mathcal{H}}{c}$", fontsize = 26)
# plt.savefig("figs/numerical_stability.pdf")
# plt.show()


# """
# Plots of Hp'/Hp and Hp''/Hp vs. x together with analytical expectations
# """
# fig, axs = plt.subplots(ncols = 2, figsize = (16, 6))
# fig.subplots_adjust(left = 0.07, right = 0.97, wspace = 0.15)

# axs[0].plot([data["x"][0], x_rm], [-1, -1], color = "yellowgreen", label = r"$r$-dominated")
# axs[0].plot([x_rm, x_mLambda], [-1/2, -1/2], color = "orange", label = r"$m$-dominated")
# axs[0].plot([x_mLambda, data["x"][-1]], [1, 1], color = "mediumorchid", label = r"$\Lambda$-dominated")
# axs[0].axvline(x_rm, 0, 1, color = "gold", label = r"$(r,m)$-equality")
# axs[0].axvline(x_acc, 0, 1, color = "coral", label = r"Onset of acceleration")
# axs[0].axvline(x_mLambda, 0, 1, color = "palevioletred", label = r"$(\Lambda,m)$-equality")
# axs[0].axvline(0, 0, 1, color = "black", linestyle = "--", label = "Today")
# axs[0].plot(data["x"], data["dHpdx"]/data["Hp"], "k", label = "Computed evolution")

# axs[1].plot([data["x"][0], x_rm], [1, 1], color = "yellowgreen", label = r"$r$-dominated")
# axs[1].plot([x_rm, x_mLambda], [1/4, 1/4], color = "orange", label = r"$m$-dominated")
# axs[1].plot([x_mLambda, data["x"][-1]], [1, 1], color = "mediumorchid", label = r"$\Lambda$-dominated")
# axs[1].axvline(x_rm, 0, 1, color = "gold", label = r"$(r,m)$-equality")
# axs[1].axvline(x_acc, 0, 1, color = "coral", label = r"Onset of acceleration")
# axs[1].axvline(x_mLambda, 0, 1, color = "palevioletred", label = r"$(\Lambda,m)$-equality")
# axs[1].axvline(0, 0, 1, color = "black", linestyle = "--", label = "Today")
# axs[1].plot(data["x"], data["ddHpddx"]/data["Hp"], "k", label = "Computed evolution")

# axs[0].legend(framealpha = 1)
# axs[0].set_xlim([x_min, x_max])
# axs[1].set_xlim([x_min, x_max])
# axs[0].set_title(r"$\frac{1}{\mathcal{H}}\frac{d\mathcal{H}}{dx}$", fontsize = 28, pad = 16)
# axs[1].set_title(r"$\frac{1}{\mathcal{H}}\frac{d^2\mathcal{H}}{dx^2}$", fontsize = 28, pad = 16)
# fig.supxlabel(r"$x=\log(a)$")
# fig.savefig("figs/H_prime_derivatives.pdf")
# plt.show()


# """
# Compute analytical expectations and shift values by constants to fit plot
# """
# x_r = np.linspace(data["x"][0], x_rm, 100)
# x_m = np.linspace(x_rm, x_mLambda, 100)
# x_Lambda = np.linspace(x_mLambda, data["x"][-1], 100)

# Hp_r = H0*np.sqrt(Omega_r0*np.exp(-2*x_r))
# Hp_m = H0*np.sqrt(Omega_m0*np.exp(-x_m))
# Hp_m += H0*(Hp_r[-1] - Hp_m[0])
# Hp_Lambda = H0*np.sqrt(Omega_Lambda0*np.exp(2*x_Lambda))
# Hp_Lambda += H0*(Hp_m[-1] - Hp_Lambda[0])

# t_r = 1/(2*H0*(1000/Mpc)) * np.exp(2*x_r)/np.sqrt(Omega_r0)
# t_m = 1/(3*H0*(1000/Mpc)) * (2*np.exp(3*x_m/2)/np.sqrt(Omega_m0) - Omega_r0**(3/2)/(2*Omega_m0**2))
# t_Lambda = 1/(H0*(1000/Mpc)*np.sqrt(Omega_Lambda0)) * (x_Lambda + 2/3 - np.log(Omega_m0/Omega_Lambda0)/3 - np.sqrt(Omega_Lambda0)*Omega_r0**(3/2)/(6*Omega_m0**2)) 

# eta_r = c/(H0*(1000/Mpc)) * np.exp(x_r)/np.sqrt(Omega_r0)
# eta_m = c/(H0*(1000/Mpc)) * (2*np.sqrt(np.exp(x_m)/Omega_m0) - np.sqrt(Omega_r0)/Omega_m0)
# eta_Lambda = - c/(H0*(1000/Mpc)) * (1/(np.sqrt(Omega_Lambda0)*np.exp(x_Lambda)) + np.sqrt(Omega_r0)/Omega_m0 - 3/(Omega_m0**(1/3)*Omega_Lambda0**(1/6)))


# """
# Plot of Hp(x) together with analytical expectations
# """
# plt.figure(figsize = (8, 6))
# plt.subplots_adjust(left = 0.09, right = 0.96, top = 0.93)
# plt.plot(x_r, Hp_r, color = "yellowgreen", linewidth = 2, linestyle = "--", label = r"$r$-dominated")
# plt.plot(x_m, Hp_m, color = "orange", linewidth = 2, linestyle = "--", label = r"$m$-dominated")
# plt.plot(x_Lambda, Hp_Lambda, color = "mediumorchid", linewidth = 2, linestyle = "--", label = r"$\Lambda$-dominated")
# plt.axvline(x_rm, 0, 1, color = "gold", label = r"$(r,m)$-equality")
# plt.axvline(x_acc, 0, 1, color = "coral", label = r"Onset of acceleration")
# plt.axvline(x_mLambda, 0, 1, color = "palevioletred", label = r"$(m,\Lambda)$-equality")
# plt.plot(data["x"], data["Hp"] * (Mpc/1000), "k", linewidth = 3, label = "Computed evolution", zorder = 0)
# plt.legend(loc = "upper right", framealpha = 1)
# plt.xlim([x_min, x_max])
# plt.yscale("log")
# plt.xlabel(r"$x=\log(a)$")
# plt.ylabel(r"$\mathcal{H}(x)$ [km/s/Mpc]")
# plt.savefig("figs/H_prime.pdf")
# plt.show()


# """
# Plot of cosmic time t(x) and conformal time eta(x)/x with analytical expectations
# """
# fig = plt.figure(layout = "constrained", figsize = (16, 12))
# subfigs = fig.subfigures(2, 1, hspace = 0.05, height_ratios = [1.5, 1])
# axs = list(subfigs[1].subplots(1, 2))
# axs.append(subfigs[0].subplots(1, 1))

# for i in range(3):
#     axs[i].plot(x_r, t_r / yr, color = "yellowgreen", linewidth = 2, linestyle = "--", label = r"$r$-dominated", zorder = 1)
#     axs[i].plot(x_m, t_m / yr, color = "orange", linewidth = 2, linestyle = "--", label = r"$m$-dominated", zorder = 1)
#     axs[i].plot(x_Lambda, t_Lambda / yr, color = "mediumorchid", linewidth = 2, linestyle = "--", label = r"$\Lambda$-dominated", zorder = 1)
#     axs[i].plot(data["x"], data["t"] / yr, "k", linewidth = 3, label = r"Cosmic time $t$", zorder = 0)

#     axs[i].plot(x_r, eta_r/c / yr, color = "yellowgreen", linewidth = 2, linestyle = "--", zorder = 1)
#     axs[i].plot(x_m, eta_m/c / yr, color = "orange", linewidth = 2, linestyle = "--", zorder = 1)
#     axs[i].plot(x_Lambda, eta_Lambda/c / yr, color = "mediumorchid", linewidth = 2, linestyle = "--", zorder = 1)

#     axs[i].axvline(x_rm, 0, 1, color = "gold", label = r"$(r,m)$-equality", zorder = 1)
#     axs[i].axvline(x_acc, 0, 1, color = "coral", label = r"Onset of acceleration", zorder = 1)
#     axs[i].axvline(x_mLambda, 0, 1, color = "palevioletred", label = r"$(m,\Lambda)$-equality", zorder = 1)

#     axs[i].plot(data["x"], data["eta"]/c / yr, "slategrey", linewidth = 3, label = r"Conformal time $\eta/c$", zorder = 0)

#     axs[i].set_yscale("log")

# axs[2].legend(ncols = 2, framealpha = 1)

# axs[2].set_xlim([x_min, x_max])
# axs[2].set_yticks([data["t"][0]/yr, 1/yr*60*60*24, 1, 1000, 1e6, 1e9, data["t"][-1]/yr], [f"{data['t'][0]:.1f} sec", "1 day", "1 yr", "1 kyr", "1 Myr", "1 Gyr", f"{data['t'][-1]/Gyr:.1f} Gyr"])
# axs[0].set_xlim([x_rm-1.5, x_rm+1.5])
# axs[0].set_ylim([10**(0.75*np.log10(t_rm / yr)), 10**(1.25*np.log10(t_rm / yr))])
# times = np.logspace(0.78*np.log10(t_rm / yr), 1.22*np.log10(t_rm / yr), 5)
# labels = [f"{time/1e3:.2f}" for time in times]
# axs[0].set_yticks(times)
# axs[0].set_yticklabels(labels)
# axs[0].set_ylabel("kyr")
# axs[1].set_xlim([x_mLambda-1.5, x_mLambda+1.5])
# axs[1].set_ylim([10**(0.9*np.log10(t_mLambda / yr)), 10**(1.1*np.log10(t_mLambda / yr))])
# times = np.logspace(0.91*np.log10(t_mLambda / yr), 1.09*np.log10(t_mLambda / yr), 5)
# labels = [f"{time/1e9:.2f}" for time in times]
# axs[1].set_yticks(times)
# axs[1].set_yticklabels(labels)
# axs[1].set_ylabel("Gyr")

# fig.supxlabel(r"$x=\log(a)$")
# fig.supylabel(r"Time")
# fig.savefig("figs/eta_and_t.pdf")
# plt.show()




# """
# Plot of density parameters
# """
# plt.figure(figsize = (8, 6))
# plt.subplots_adjust(left = 0.09, right = 0.96, top = 0.93)

# plt.plot(data["x"], data["Omega_gamma"] + data["Omega_nu"], color = "mediumseagreen", label = r"$\Omega_{r}$", zorder = 1)
# plt.plot(data["x"], data["Omega_gamma"], color = "yellowgreen", linestyle = "--", label = r"$\Omega_{\gamma}$", zorder = 0)
# plt.plot(data["x"], data["Omega_nu"], color = "skyblue", linestyle = "--", label = r"$\Omega_{\nu}$", zorder = 0)
# plt.plot(data["x"], data["Omega_b"] + data["Omega_CDM"], color = "orange", label = r"$\Omega_{m}$", zorder = 1)
# plt.plot(data["x"], data["Omega_b"], color = "#ffd726", linestyle = "--", label = r"$\Omega_{b}$", zorder = 0)
# plt.plot(data["x"], data["Omega_CDM"], color = "palevioletred", linestyle = "--", label = r"$\Omega_{\small\textrm{CDM}}$", zorder = 0)
# plt.plot(data["x"], data["Omega_Lambda"], color = "mediumorchid", label = r"$\Omega_{\Lambda}$")

# plt.legend(loc = "center left", framealpha = 1)
# plt.xlim([x_min, x_max])
# plt.xlabel(r"$x=\log(a)$")
# plt.ylabel(r"Density parameters $\Omega_i$")
# plt.savefig("figs/density_parameters.pdf")
# plt.show()
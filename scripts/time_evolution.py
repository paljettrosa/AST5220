import numpy as np
import matplotlib.pyplot as plt
plt.rcParams.update({'text.usetex': True, 'font.size': 18, 'font.family': 'serif', 'font.serif': 'Computer Modern Sans Serif', 'font.weight': 100, 'mathtext.fontset': 'cm', 'xtick.labelsize': 16, 'ytick.labelsize': 16})

#TODO remove this
c = 299792458
hbar = 1.05457182e-34
k_B = 1.380649e-23
T_CMB0 = 2.725
m3_to_cm3 = 1e6
zeta_3 = 1.2020569031595942853

n_nu0 = 3*zeta_3/(2*np.pi**2) * 4/11 * (T_CMB0*k_B/(hbar*c))**3
print(n_nu0 / m3_to_cm3)


"""
Calculate analytical expressions
"""
c = 299792458               # [m/s]
Mpc = 3.08567758e22         # [m]
yr = 3.1536e7               # [s]
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
t_rm = Omega_r0**(3/2) / (2*H0*Omega_m0**2) * (Mpc/1000)
print("Radiation-matter equality:")
print(f"x: {x_rm:.3f}")
print(f"z: {z_rm:.3f}")
print(f"t: {t_rm/Myr:.3f} Myr\n")

# Matter-dark energy equality
x_mLambda = np.log(Omega_m0/Omega_Lambda0) / 3
z_mLambda = (Omega_Lambda0/Omega_m0)**(1/3) - 1
t_mLambda = 1/(3*H0) * (2/np.sqrt(Omega_Lambda0) - Omega_r0**(3/2) / (2*Omega_m0**2)) * (Mpc/1000)
print("Matter-dark energy equality:")
print(f"x: {x_mLambda:.3f}")
print(f"z: {z_mLambda:.3f}")
print(f"t: {t_mLambda/Gyr:.3f} Gyr\n")

# Onset of acceleration
x_acc = np.log(Omega_m0/(2*Omega_Lambda0)) / 3
z_acc = (2*Omega_Lambda0/Omega_m0)**(1/3) - 1
t_acc = 1/(3*H0) * (np.sqrt(2/Omega_Lambda0) - Omega_r0**(3/2) / (2*Omega_m0**2)) * (Mpc/1000)
print("Onset of acceleration:")
print(f"x: {x_acc:.3f}")
print(f"z: {z_acc:.3f}")
print(f"t: {t_acc/Gyr:.3f} Gyr\n")


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

TODO maybe compute eta*Hp/c and make plot. add analytical expectations. have all four tests in one figure?
"""
# plt.figure(figsize = (8, 6))
# plt.subplots_adjust(left = 0.16, right = 0.95, top = 0.93)
# plt.plot(data["x"], data["detadx"]*data["Hp"]/c, "k")
# plt.plot([data["x"][0], data["x"][-1]], [1, 1], color = "#de82b4", linestyle = "--")
# plt.xlabel(r"$x=\log(a)$")
# plt.ylabel(r"$\frac{\eta'\mathcal{H}}{c}$", fontsize = 24)
# plt.savefig("figs/numerical_stability.pdf")
# plt.show()


"""
Plots of Hp'/Hp and Hp''/Hp vs. x together with analytical expectations
"""
# fig, axs = plt.subplots(ncols = 2, figsize = (16, 6))
# fig.subplots_adjust(left = 0.07, right = 0.97, wspace = 0.15)

# axs[0].plot([data["x"][0], x_rm], [-1, -1], color = "yellowgreen", label = r"$r$-dominated")
# axs[0].plot([x_rm, x_rm], [-1, -1/2], color = "gold", linestyle = "--")
# axs[0].plot([x_rm, x_mLambda], [-1/2, -1/2], color = "orange", label = r"$m$-dominated")
# axs[0].plot([x_mLambda, x_mLambda], [-1/2, 1], color = "palevioletred", linestyle = "--")
# axs[0].plot([x_mLambda, data["x"][-1]], [1, 1], color = "mediumorchid", label = r"$\Lambda$-dominated")
# axs[0].plot(data["x"], data["dHpdx"]/data["Hp"], "k", label = "Computed evolution")

# axs[1].plot([data["x"][0], x_rm], [1, 1], color = "yellowgreen", label = r"$r$-dominated")
# axs[1].plot([x_rm, x_rm], [1, 1/4], color = "gold", linestyle = "--")
# axs[1].plot([x_rm, x_mLambda], [1/4, 1/4], color = "orange", label = r"$m$-dominated")
# axs[1].plot([x_mLambda, x_mLambda], [1/4, 1], color = "palevioletred", linestyle = "--")
# axs[1].plot([x_mLambda, data["x"][-1]], [1, 1], color = "mediumorchid", label = r"$\Lambda$-dominated")
# axs[1].plot(data["x"], data["ddHpddx"]/data["Hp"], "k", label = "Computed evolution")

# axs[0].legend()
# axs[0].set_title(r"$\frac{1}{\mathcal{H}}\frac{d\mathcal{H}}{dx}$", fontsize = 28, pad = 16)
# axs[1].set_title(r"$\frac{1}{\mathcal{H}}\frac{d^2\mathcal{H}}{dx^2}$", fontsize = 28, pad = 16)
# fig.supxlabel(r"$x=\log(a)$")
# fig.savefig("figs/H_prime_derivatives.pdf")
# plt.show()





"""
Compute analytical expectations for the scale factor, and shift values by constants to fit plot

#TODO maybe do differently. spline x(t)?
#TODO express things in terms of equality times instead?
"""
# t_r = np.linspace(data["t"][0], t_rm, 100)
# t_m = np.linspace(t_rm, t_mLambda, 100)
# t_Lambda = np.linspace(t_mLambda, data["t"][-1], 100)

# x_r = np.log(np.sqrt(t_r))
# x_m = np.log((t_m)**(2/3)) 
# x_Lambda = H0 * (1000/Mpc) * np.sqrt(Omega_Lambda0) * t_Lambda

# x_r -= x_r[0] - data["x"][0]
# x_m -= x_m[0] - x_rm
# x_Lambda -= x_Lambda[-1] - data["x"][-1]

x_r = np.linspace(data["x"][0], x_rm, 100)
x_m = np.linspace(x_rm, x_mLambda, 100)
x_Lambda = np.linspace(x_mLambda, data["x"][-1], 100)

t_r = 1/(2*H0*(1000/Mpc)) * np.exp(2*x_r)/np.sqrt(Omega_r0)
t_m = 1/(3*H0*(1000/Mpc)) * (2*np.exp(3*x_m/2)/np.sqrt(Omega_m0) - Omega_r0**(3/2)/(2*Omega_m0**2))
t_Lambda = 1/(H0*(1000/Mpc)*np.sqrt(Omega_Lambda0)) * (x_Lambda + 2/3 - np.log(Omega_m0/Omega_Lambda0)/3 - np.sqrt(Omega_Lambda0)*Omega_r0**(3/2)/(6*Omega_m0**2)) 

eta_r = c/(H0*(1000/Mpc)) * np.exp(x_r)/np.sqrt(Omega_r0)
eta_m = c/(H0*(1000/Mpc)) * (2*np.sqrt(np.exp(x_m)/Omega_m0) - np.sqrt(Omega_r0)/Omega_m0)
eta_Lambda = - c/(H0*(1000/Mpc)) * (1/(np.sqrt(Omega_Lambda0)*np.exp(x_Lambda)) + np.sqrt(Omega_r0)/Omega_m0 - 3/(Omega_m0**(1/3)*Omega_Lambda0**(1/6)))

"""
Plot of cosmic time t(x) and conformal time eta(x)/x with analytical expectations

#TODO visualize in another way?
#TODO use another time unit?
#TODO analytical expectations for eta as well?
#TODO argue that it is inaccurate in the matter dominated era since radiation is not negligible in the beginning and dark energy not in the end?
#TODO insert age of the universe today
"""
# fig, axs = plt.subplots(ncols = 2, figsize = (16, 6))
# fig.subplots_adjust(left = 0.07, right = 0.97, wspace = 0.15)

# axs[0].plot(x_r, t_r / yr, color = "yellowgreen", linestyle = "--", label = r"$r$-dominated", zorder = 1)
# axs[0].plot([data["x"][0], data["x"][-1]], [t_rm / yr, t_rm / yr], color = "gold", label = r"$r-m$-equality", zorder = 1)
# axs[0].plot(x_m, t_m / yr, color = "orange", linestyle = "--", label = r"$m$-dominated", zorder = 1)
# axs[0].plot([data["x"][0], data["x"][-1]], [t_mLambda / yr, t_mLambda / yr], color = "palevioletred", label = r"$m-\Lambda$-equality", zorder = 1)
# axs[0].plot(x_Lambda, t_Lambda / yr, color = "mediumorchid", linestyle = "--", label = r"$\Lambda$-dominated", zorder = 1)
# axs[0].plot(data["x"], data["t"] / yr, "k", label = "Computed evolution", zorder = 0)

# axs[1].plot(data["x"], data["eta"]/c / yr, "k")

# axs[0].legend(loc = "lower right")
# axs[0].set_yscale("log")
# axs[1].set_yscale("log")
# fig.supxlabel(r"$x=\log(a)$")
# fig.supylabel(r"Time [yr]")
# axs[0].set_title(r"Cosmic time $t$")
# axs[1].set_title(r"Conformal time $\eta/c$")
# fig.savefig("figs/eta_and_t.pdf")
# plt.show()

plt.figure(figsize = (8, 6))
plt.subplots_adjust(left = 0.17, right = 0.96, top = 0.93)

plt.plot(x_r, t_r / yr, color = "yellowgreen", linewidth = 2, linestyle = "--", label = r"$r$-dominated", zorder = 1)
plt.plot(x_m, t_m / yr, color = "orange", linewidth = 2, linestyle = "--", label = r"$m$-dominated", zorder = 1)
plt.plot(x_Lambda, t_Lambda / yr, color = "mediumorchid", linewidth = 2, linestyle = "--", label = r"$\Lambda$-dominated", zorder = 1)
plt.plot(data["x"], data["t"] / yr, "k", linewidth = 3, label = r"Cosmic time $t$", zorder = 0)

plt.plot(x_r, eta_r/c / yr, color = "yellowgreen", linewidth = 2, linestyle = "--", zorder = 1)
plt.plot(x_m, eta_m/c / yr, color = "orange", linewidth = 2, linestyle = "--", zorder = 1)
plt.plot(x_Lambda, eta_Lambda/c / yr, color = "mediumorchid", linewidth = 2, linestyle = "--", zorder = 1)
plt.plot(data["x"], data["eta"]/c / yr, "slategrey", linewidth = 3, label = r"Conformal time $\eta/c$", zorder = 0)

plt.axvline(x_rm, 0, 1, color = "gold", label = r"$(r,m)$-equality", zorder = 1)
plt.axvline(x_mLambda, 0, 1, color = "palevioletred", label = r"$(m,\Lambda)$-equality", zorder = 1)

t_0 = 13.7e9 #TODO
plt.legend(framealpha = 1)
plt.yscale("log")
plt.xlabel(r"$x=\log(a)$")
plt.ylabel(r"Time")
plt.yticks([data["t"][0]/yr, 1/yr*60*60*24, 1, 1000, 1e6, 1e9, data["t"][-1]/yr], [f"{data['t'][0]:.1f} sec", "1 day", "1 yr", "1000 yr", "1 Myr", "1 Gyr", f"{data['t'][-1]/Gyr:.1f} Gyr"])
plt.savefig("figs/eta_and_t.pdf")
plt.show()



"""
Plot of x vs. cosmic time t

#TODO maybe not necessary?
"""
# plt.plot(t_r / yr, x_r, color = "yellowgreen", linewidth = 5)
# plt.plot([t_rm / yr, t_rm / yr], [data["x"][0], data["x"][-1]], color = "gold")
# plt.plot(t_m / yr, x_m, color = "orange", linewidth = 5)
# plt.plot([t_mLambda / yr, t_mLambda / yr], [data["x"][0], data["x"][-1]], color = "palevioletred")
# plt.plot(t_Lambda / yr, x_Lambda, color = "mediumorchid", linewidth = 5)
# plt.plot(data["t"] / yr, data["x"], "k")
# plt.xscale("log")
# plt.show()



"""
Plot of density parameters
"""
plt.figure(figsize = (8, 6))
plt.subplots_adjust(left = 0.09, right = 0.96, top = 0.93)

plt.plot(data["x"], data["Omega_gamma"] + data["Omega_nu"], color = "yellowgreen", label = r"$\Omega_{r} = \Omega_{\gamma} + \Omega_{\nu}$")
plt.plot(data["x"], data["Omega_b"] + data["Omega_CDM"], color = "orange", label = r"$\Omega_{m} = \Omega_{b} + \Omega_{\small\textrm{CDM}}$")
plt.plot(data["x"], data["Omega_Lambda"], color = "mediumorchid", label = r"$\Omega_{\Lambda}$")

# plt.plot(data["x"], data["Omega_gamma"] + data["Omega_nu"], color = "mediumseagreen", label = r"$\Omega_{r}$", zorder = 1)
# plt.plot(data["x"], data["Omega_gamma"], color = "yellowgreen", linestyle = "--", label = r"$\Omega_{\gamma}$", zorder = 0)
# plt.plot(data["x"], data["Omega_nu"], color = "skyblue", linestyle = "--", label = r"$\Omega_{\nu}$", zorder = 0)
# plt.plot(data["x"], data["Omega_b"] + data["Omega_CDM"], color = "orange", label = r"$\Omega_{m}$", zorder = 1)
# plt.plot(data["x"], data["Omega_b"], color = "#ffd726", linestyle = "--", label = r"$\Omega_{b}$", zorder = 0)
# plt.plot(data["x"], data["Omega_CDM"], color = "palevioletred", linestyle = "--", label = r"$\Omega_{\small\textrm{CDM}}$", zorder = 0)
# plt.plot(data["x"], data["Omega_Lambda"], color = "mediumpurple", label = r"$\Omega_{\Lambda}$")

plt.legend()
plt.xlabel(r"$x=\log(a)$")
plt.ylabel(r"Density parameters $\Omega_i$")
plt.savefig("figs/density_parameters.pdf")
plt.show()

#TODO integrate from -21 to 6 instead, or something like that
#TODO use points instead of lines in the eta derivative plot
#TODO maybe change font to astro font, and figsize to ish nine, fontsize to text fontsize
#TODO maybe make some subfigures showing deviations from analytical expressions for the eta/t plot. two tall subplots next to the normal plot (zoomed in on r-m and m-Lambda equalities), or own figure? use linear time axis in that case?
#TODO fix supernova plots. add distribution function to histograms (print standard deviations). maybe only plot H0?
#TODO calculate all times asked for, use spline for accurate value of time today
#TODO write intro and theory done
#TODO write implementation and testing, move some theory to this section?
#TODO maybe do neutrino thing?

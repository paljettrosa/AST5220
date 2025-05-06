import numpy as np
import matplotlib.pyplot as plt
plt.rcParams.update({'text.usetex': True, 'font.size': 18, 'font.family': 'serif', 'font.serif': 'Computer Modern Sans Serif', 'font.weight': 100, 'mathtext.fontset': 'cm', 'xtick.labelsize': 16, 'ytick.labelsize': 16})

"""
####################################
Milestone I: Background Cosmology
####################################
"""
c = 299792458               # [m/s]
Mpc = 3.08567758e22         # [m]

x_min = -20
x_max = 5
n_pts = (x_max - x_min)*100 + 1
data = {}
keys = ["x", "eta", "Hp", "dHpdx", "ddHpddx", "Omega_b", "Omega_CDM", "Omega_Lambda", "Omega_gamma", "Omega_nu", "Omega_k"]
for j, key in enumerate(keys):
    data[key] = np.zeros(n_pts)

with open("results/toy/cosmology.txt", "r") as infile:
    for i, line in enumerate(infile):
        values = line.split()
        for j, key in enumerate(keys):
            data[key][i] = float(values[j])

"""
Hubble parameter and conformal time
"""
fig, axs = plt.subplots(ncols = 2, layout = "constrained", figsize = (16, 6))

axs[0].plot(data["x"], data["Hp"]*Mpc/1e5, color = "cornflowerblue", linewidth = 2)
axs[1].plot(data["x"], data["eta"]/Mpc, color = "palevioletred", linewidth = 2)

axs[0].set_yscale("log")
axs[1].set_yscale("log")
axs[0].set_ylim(1e-1, 1e3)
axs[1].set_ylim(1, 3e4)
axs[0].set_xlim(-12, 0)
axs[1].set_xlim(-12, 0)
axs[0].set_ylabel(r"$\mathcal{H}$ [$100\,$km/s/Mpc]")
axs[1].set_ylabel(r"$\eta$ [Mpc]")

fig.supxlabel(r"$x=\log(a)$")
fig.savefig("figs/toy/Hp_and_eta.pdf")
plt.show()


"""
Fraction eta*Hp/c
"""
plt.figure(figsize = (8, 6))
plt.subplots_adjust(left = 0.16, right = 0.95, top = 0.93)
plt.plot(data["x"], data["eta"]*data["Hp"]/c, color = "mediumpurple", linewidth = 2)
plt.xlim(-15, 0)
plt.ylim(0.8, 3)
plt.xlabel(r"$x=\log(a)$")
plt.ylabel(r"$\frac{\eta\mathcal{H}}{c}$", fontsize = 26)
plt.savefig("figs/toy/eta_Hp_over_c.pdf")
plt.show()


"""
Density parameters
"""
plt.figure(figsize = (8, 6))
plt.subplots_adjust(left = 0.09, right = 0.96, top = 0.93)

plt.plot(data["x"], data["Omega_gamma"] + data["Omega_nu"], color = "mediumseagreen", label = r"$\Omega_{r}=\Omega_\gamma+\Omega_\nu$")
plt.plot(data["x"], data["Omega_b"] + data["Omega_CDM"], color = "orange", label = r"$\Omega_{m}=\Omega_b+\Omega_{\small\textrm{CDM}}$")
plt.plot(data["x"], data["Omega_Lambda"], color = "mediumorchid", label = r"$\Omega_{\Lambda}$")

plt.legend(loc = "center left", framealpha = 1)
plt.xlim(x_min, x_max)
plt.ylim(0, 1.1)
plt.xlabel(r"$x=\log(a)$")
plt.ylabel(r"Density parameters $\Omega_i$")
plt.savefig("figs/toy/Omega_i.pdf")
plt.show()



"""
####################################
Milestone II: Recombination History
####################################
"""
x_min = -12
x_max = 0
n_pts = (x_max - x_min)*10000 + 1

data = {}
data_reion = {}
data_Hereion = {}
keys = ["x", "Xe", "ne", "tau", "dtaudx", "ddtauddx", "g_tilde", "dgdx_tilde", "ddgddx_tilde"]
for key in keys:
    data[key] = np.zeros(n_pts)
    data_reion[key] = np.zeros(n_pts)
    data_Hereion[key] = np.zeros(n_pts)

with open("results/toy/recombination.txt", "r") as infile:
    for i, line in enumerate(infile):
        values = line.split()
        for j, key in enumerate(keys):
            data[key][i] = float(values[j])

with open("results/toy/recombination_reion.txt", "r") as infile:
    for i, line in enumerate(infile):
        values = line.split()
        for j, key in enumerate(keys):
            data_reion[key][i] = float(values[j])

with open("results/toy/recombination_Hereion.txt", "r") as infile:
    for i, line in enumerate(infile):
        values = line.split()
        for j, key in enumerate(keys):
            data_Hereion[key][i] = float(values[j])

"""
Free electron fraction
"""
fig, axs = plt.subplots(ncols = 2, layout = "constrained", figsize = (16, 6))

axs[0].plot(data["x"], data["Xe"], color = "#b55087", label = r"Without reionization", linewidth = 2)
axs[0].plot(data_reion["x"], data_reion["Xe"], color = "#de82b4", linestyle = "--", label = r"With reionization", linewidth = 2)
axs[1].plot(data_Hereion["x"], data_Hereion["Xe"], color = "#b55087", linewidth = 2)

axs[0].legend(loc = "upper right", framealpha = 1)
axs[0].set_yscale("log")
axs[0].set_ylim(1e-4, 1e1)
axs[0].set_xlim(x_min, x_max)
axs[1].set_xlim(x_min, x_max)
axs[0].set_title(r"Excluding Helium")
axs[1].set_title(r"Including Helium and reionization")

fig.supxlabel(r"$x=\log(a)$")
fig.supylabel(r"Free electron fraction $X_e$")
fig.savefig("figs/toy/Xe.pdf")
plt.show()


"""
Optical depth
"""
fig, axs = plt.subplots(ncols = 2, layout = "constrained", figsize = (16, 6))

axs[0].plot(data["x"], data["tau"], color = "cornflowerblue", linewidth = 2)
axs[0].plot(data["x"], -data["dtaudx"], color = "palevioletred", linewidth = 2)
axs[0].plot(data["x"], data["ddtauddx"], color = "yellowgreen", linewidth = 2)
axs[1].plot(data_reion["x"], data_reion["tau"], color = "cornflowerblue", label = r"$\tau(x)$", linewidth = 2)
axs[1].plot(data_reion["x"], -data_reion["dtaudx"], color = "palevioletred", label = r"-$\tau'(x)$", linewidth = 2)
axs[1].plot(data_reion["x"], data_reion["ddtauddx"], color = "yellowgreen", label = r"$\tau''(x)$", linewidth = 2)

axs[1].legend(loc = "upper right", framealpha = 1)
axs[0].set_yscale("log")
axs[1].set_yscale("log")
axs[0].set_ylim(1e-8, 1e8)
axs[1].set_ylim(1e-8, 1e8)
axs[0].set_xlim(x_min, x_max)
axs[1].set_xlim(x_min, x_max)
axs[0].set_title(r"Without reionization")
axs[1].set_title(r"With reionization")

fig.supxlabel(r"$x=\log(a)$")
fig.supylabel(r"Optical depth and derivatives")
fig.savefig("figs/toy/tau.pdf")
plt.show()


"""
Visibility function
"""
fig, axs = plt.subplots(nrows = 3, layout = "constrained", figsize = (8, 14))

quantities = [data["g_tilde"], data["dgdx_tilde"], data["ddgddx_tilde"], data_reion["g_tilde"], data_reion["dgdx_tilde"], data_reion["ddgddx_tilde"]]
titles = [r"$\tilde{g}(x)$", r"$\tilde{g}'(x)$", r"$\tilde{g}''(x)$"]
for i in range(3):
    axs[i].plot(data["x"], quantities[i], "crimson", label = "Without reionization")
    axs[i].plot(data_reion["x"], quantities[i+3], "pink", label = "With reionization")
    axs[i].set_title(titles[i])
    axs[i].set_xlim(x_min , x_max)
axs[0].legend(loc = "upper right")

fig.supxlabel(r"$x=\log(a)$")
fig.supylabel(r"Visibility function and derivatives")
fig.savefig("figs/toy/g_tilde.pdf")
plt.show()


"""
####################################
Milestone III: Perturbations
####################################
"""
x_min = -18
x_max = 0
n_pts = (x_max - x_min)*10000 + 1

data = {}
k_vals = [0.001, 0.01, 0.1]
keys = ["x", "Phi", "Psi", "Pi", "delta_CDM", "delta_b", "v_CDM", "v_b", "Theta_0", "Theta_1", "Theta_2"]
for k in k_vals:
    data[k] = {}
    for key in keys:
        data[k][key] = np.zeros(n_pts)

    with open(f"results/toy/perturbations_k{k}.txt", "r") as infile:
        for i, line in enumerate(infile):
            values = line.split()
            for j, key in enumerate(keys):
                data[k][key][i] = float(values[j])

colors = ["cornflowerblue", "palevioletred", "yellowgreen"]

"""
Matter densities and velocities
"""
fig, axs = plt.subplots(ncols = 2, layout = "constrained", figsize = (16, 6))

for i, k in enumerate(k_vals):
    axs[0].plot(data[k]["x"], abs(data[k]["delta_CDM"]), color = colors[i], label = r"$k=\,$" + f"{k}" + r"$\,\textrm{Mpc}^{-1}$")
    axs[0].plot(data[k]["x"], abs(data[k]["delta_b"]), color = colors[i], linestyle = "--")
    axs[1].plot(data[k]["x"], abs(data[k]["v_CDM"]), color = colors[i])
    axs[1].plot(data[k]["x"], abs(data[k]["v_b"]), color = colors[i], linestyle = "--") 

axs[0].legend(loc = "upper left", framealpha = 1)
axs[0].set_yscale("log")
axs[1].set_yscale("log")
axs[0].set_ylim(1e-1, 1e5)
axs[0].set_title(r"Densities $\delta_{\small\textrm{CDM}}$ (solid) and $\delta_b$ (dashed)")
axs[1].set_title(r"Velocities $v_{\small\textrm{CDM}}$ (solid) and $v_b$ (dashed)")

fig.supxlabel(r"$x=\log(a)$")
fig.savefig("figs/toy/delta_and_v.pdf")
plt.show()


"""
Photon temperature monopole and dipole
"""
fig, axs = plt.subplots(ncols = 2, layout = "constrained", figsize = (16, 6))

for i, k in enumerate(k_vals):
    axs[0].plot(data[k]["x"], data[k]["Theta_0"], color = colors[i], label = r"$k=\,$" + f"{k}" + r"$\,\textrm{Mpc}^{-1}$")
    axs[1].plot(data[k]["x"], data[k]["Theta_1"], color = colors[i]) 

axs[0].legend(loc = "lower left", framealpha = 1)
axs[0].set_title(r"Monopole $\Theta_0$")
axs[1].set_title(r"Dipole $\Theta_1$")

fig.supxlabel(r"$x=\log(a)$")
fig.savefig("figs/toy/Thetas.pdf")
plt.show()


"""
Gravitational potential
"""
plt.figure(figsize = (8, 6))
plt.subplots_adjust(left = 0.09, right = 0.96, top = 0.93)

for i, k in enumerate(k_vals):
    plt.plot(data[k]["x"], data[k]["Phi"], color = colors[i], label = r"$k=\,$" + f"{k}" + r"$\,\textrm{Mpc}^{-1}$")

plt.legend(loc = "lower left", framealpha = 1)
plt.xlabel(r"$x=\log(a)$")
plt.ylabel(r"Gravitational potential $\Phi$")
plt.savefig("figs/toy/Phi.pdf")
plt.show()



"""
####################################
Milestone IV: Power Spectra
####################################
"""
ell_max = 2000
n_ells  = ell_max-1

k_min   = 0.00005
k_max   = 0.3
n_pts_k = int(np.log(k_max) - np.log(k_min))*10000 + 1

n_pts   = [n_ells, n_pts_k]

data      = {}
filenames = ["C_TT", "P_k"]
keys      = [["ell", "C_TT"], ["k", "P_k"]]
for i, filename in enumerate(filenames):
    data[filename] = {}
    for key in keys[i]:
        data[filename][key] = np.zeros(n_pts[i])

    with open(f"results/toy/{filename}.txt", "r") as infile:
        for j, line in enumerate(infile):
            values = line.split()
            for k, key in enumerate(keys[i]):
                data[filename][key][j] = float(values[k])


fig, axs = plt.subplots(ncols = 2, layout = "constrained", figsize = (16, 6))

axs[0].plot(data["C_TT"]["ell"], data["C_TT"]["C_TT"], color = "#d767ad", linewidth = 2)
axs[1].plot(data["P_k"]["k"], data["P_k"]["P_k"], color = "#ffab0f", linewidth = 2)

axs[0].set_xscale("log")
axs[1].set_xscale("log")
axs[0].set_yscale("log")
axs[1].set_yscale("log")
axs[0].set_xlabel(r"Multipole $\ell$")
axs[1].set_xlabel(r"Wavenumber $k$ [$h/\textrm{Mpc}$]")
axs[0].set_ylabel(r"$\ell(\ell+1)C_\ell/2\pi$ [$\mu\textrm{K}^2$]")
axs[1].set_ylabel(r"$P(k)$ [$(\textrm{Mpc}/h)^3$]")
axs[0].set_title(r"Photon temperature power spectrum")
axs[1].set_title(r"Matter power spectrum")

fig.savefig("figs/toy/C_ell_and_P_k.pdf")
plt.show()
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.legend_handler import HandlerTuple
import matplotlib.patheffects as pe
plt.rcParams.update({'text.usetex': True, 'font.size': 18, 'font.family': 'serif', 'font.serif': 'Computer Modern Sans Serif', 'font.weight': 100, 'mathtext.fontset': 'cm', 'xtick.labelsize': 16, 'ytick.labelsize': 16})

"""
Define important times
"""
# Density parameters today
Omega_m0 = 0.05 + 0.267
Omega_r0 = 5.50896e-05 + 3.81093e-05
Omega_Lambda0 = 0.682907

# Radiation-matter equality
x_rm = np.log(Omega_r0/Omega_m0) 
# Matter-dark energy equality
x_mLambda = np.log(Omega_m0/Omega_Lambda0) / 3      


"""
Read data from files
"""
x_min = -18
x_max = 0
n_pts = (x_max - x_min)*10000 + 1

data = {}
k_vals = [1.0, 0.1, 0.01, 0.001]
keys = ["x", "Phi", "Psi", "Pi", "delta_CDM", "delta_b", "v_CDM", "v_b", "Theta_0", "Theta_1", "Theta_2", "Theta_P0", "Theta_P1", "Theta_P2", "Nu_0", "Nu_1", "Nu_2"]
for k in k_vals:
    data[k] = {}
    for key in keys:
        data[k][key] = np.zeros(n_pts)

    with open(f"results/perturbations_k{k}.txt", "r") as infile:
        for i, line in enumerate(infile):
            values = line.split()
            for j, key in enumerate(keys):
                data[k][key][i] = float(values[j])

colors_grav  = ["#f6bad0", "#eb6395", "#d91b62", "#901241"]
colors_grav  = ["#f6bad0", "#eb6395", "#c90e53", "#820131"]
colors_CDM   = ["#dfc1d6", "#ad85b0", "#825685", "#523654"]
colors_CDM   = ["#dfb1e3", "#b87bbd", "#8a4a8f", "#4f3052"]
colors_b     = ["#f7cb79", "#f7a305", "#c37e00", "#825400"]
colors_b     = ["#fad184", "#f09d02", "#a66c02", "#825400"]
colors_gamma = ["#e6e995", "#bec427", "#90951e", "#5b5f13"]
colors_gamma = ["#e1e38f", "#a8ad17", "#7a8018", "#51540e"]
colors_nu    = ["#b2c8e5", "#759cd0", "#3d71b7", "#294b7a"]
colors_nu    = ["#b2c8e5", "#6791c7", "#2f5e9c", "#1e3c66"]
colors_pol   = ["#fdbaa2", "#fc7341", "#bf3604", "#902702"]
colors_pol   = ["#fdbaa2", "#fa6934", "#bf3604", "#902702"]

x_enter = [-13.04, -10.71, -8.08, -4.13] 


"""
Gravitational potentials Phi and Phi + Psi
"""
fig, axs = plt.subplots(ncols = 2, layout = "constrained", figsize = (16, 6))

for i, k in enumerate(k_vals):
    axs[0].plot(data[k]["x"], data[k]["Phi"], colors_grav[i], label = r"$k=\,$" + f"{k}", linewidth = 1)
    axs[1].plot(data[k]["x"], data[k]["Phi"] + data[k]["Psi"], colors_grav[i], linewidth = 1)
    axs[0].axvline(x_enter[i], color = colors_grav[i], linestyle = "--", linewidth = 1)
    axs[1].axvline(x_enter[i], color = colors_grav[i], linestyle = "--", linewidth = 1)
for i in range(2):
    axs[i].axvline(x_rm, color = "gold", label = r"$(r,m)$-eq." if i == 0 else None, linewidth = 1, path_effects = [pe.Stroke(linewidth = 2, foreground = "black"), pe.Normal()])
    axs[i].axvline(x_mLambda, color = "palevioletred", label = r"$(m,\Lambda)$-eq." if i == 0 else None, linewidth = 1, path_effects = [pe.Stroke(linewidth = 2, foreground = "black"), pe.Normal()])

axs[0].legend(loc = "lower left", framealpha = 1)
axs[0].set_xlim(x_min, x_max)
axs[1].set_xlim(x_min, x_max)
axs[0].set_title(r"Gravitational potential $\Phi$")
axs[1].set_title(r"Anisotropic stress $\Phi+\Psi$")
fig.supxlabel(r"$x=\log(a)$")
fig.savefig("figs/Phi_and_Psi.pdf")
plt.show()


"""
CDM and baryon density and velocity perturbations
"""
fig, axs = plt.subplots(ncols = 2, layout = "constrained", figsize = (16, 6))

lines = []
labels = []
for i, k in enumerate(k_vals):
    l1 = axs[0].plot(data[k]["x"], data[k]["delta_b"], colors_b[i], linewidth = 1)[0]
    l2 = axs[0].plot(data[k]["x"], data[k]["delta_CDM"], colors_CDM[i], linewidth = 1)[0]
    lines.append((l1, l2))
    labels.append(r"$k=\,$" + f"{k}")
    axs[1].plot(data[k]["x"], data[k]["v_b"], colors_b[i], linewidth = 1)
    axs[1].plot(data[k]["x"], data[k]["v_CDM"], colors_CDM[i], linewidth = 1)
    for j in range(2):
        axs[j].axvline(x_enter[i], color = colors_b[i], linewidth = 1)
        axs[j].axvline(x_enter[i], color = colors_CDM[i], linestyle = "--", linewidth = 1)

lines.append(axs[0].axvline(x_rm, color = "gold", linewidth = 1, path_effects = [pe.Stroke(linewidth = 2, foreground = "black"), pe.Normal()]))
labels.append(r"$(r,m)$-eq.")
axs[1].axvline(x_rm, color = "gold", linewidth = 1, path_effects = [pe.Stroke(linewidth = 2, foreground = "black"), pe.Normal()])
lines.append(axs[0].axvline(x_mLambda, color = "palevioletred", linewidth = 1, path_effects = [pe.Stroke(linewidth = 2, foreground = "black"), pe.Normal()]))
labels.append(r"$(\Lambda,m)$-eq.")
axs[1].axvline(x_mLambda, color = "palevioletred", linewidth = 1, path_effects = [pe.Stroke(linewidth = 2, foreground = "black"), pe.Normal()])

axs[0].legend(lines, labels, handler_map = {tuple: HandlerTuple(ndivide = None)}, loc = "upper left", framealpha = 1)
axs[0].set_xlim(x_min, x_max)
axs[1].set_xlim(x_min, x_max)
axs[0].set_ylim(bottom = -8)
axs[0].set_yscale("asinh")
axs[1].set_yscale("asinh")
axs[0].set_yticks([1e4, 1e3, 1e2, 1e1, 2, 0, -2], labels = [r"$10^4$", r"$10^3$", "100", "10", "2", "0", "-2"])
axs[1].set_yticks([20, 10, 5, 2, 0.5, 0, -0.5, -2], labels = ["20", "10", "5", "2", "0.5", "0", "-0.5", "-2"])
axs[0].set_title(r"Densities $\delta_{\small\textrm{CDM}}$ (purple) and $\delta_b$ (orange)")
axs[1].set_title(r"Velocities $v_{\small\textrm{CDM}}$ (purple) and $v_b$ (orange)")
fig.supxlabel(r"$x=\log(a)$")
fig.savefig("figs/CDM_and_baryons.pdf")
plt.show()


"""
Photon density and velocity perturbations
"""
fig, axs = plt.subplots(ncols = 2, layout = "constrained", figsize = (16, 6))

for i, k in enumerate(k_vals):
    axs[0].plot(data[k]["x"], 4*data[k]["Theta_0"], colors_gamma[i], label = r"$k=\,$" + f"{k}", linewidth = 1)
    axs[1].plot(data[k]["x"], -3*data[k]["Theta_1"], colors_gamma[i], linewidth = 1)
    axs[0].axvline(x_enter[i], color = colors_gamma[i], linestyle = "--", linewidth = 1)
    axs[1].axvline(x_enter[i], color = colors_gamma[i], linestyle = "--", linewidth = 1)
for i in range(2):
    axs[i].axvline(x_rm, color = "gold", label = r"$(r,m)$-eq." if i == 0 else None, linewidth = 1, path_effects = [pe.Stroke(linewidth = 2, foreground = "black"), pe.Normal()])
    axs[i].axvline(x_mLambda, color = "palevioletred", label = r"$(m,\Lambda)$-eq." if i == 0 else None, linewidth = 1, path_effects = [pe.Stroke(linewidth = 2, foreground = "black"), pe.Normal()])

axs[0].legend(loc = "lower left", framealpha = 1)
axs[0].set_xlim(x_min, x_max)
axs[1].set_xlim(x_min, x_max)
axs[0].set_title(r"Density $\delta_\gamma=4\Theta_0$")
axs[1].set_title(r"Velocity $v_\gamma=-3\Theta_1$")
fig.supxlabel(r"$x=\log(a)$")
fig.savefig("figs/photons.pdf")
plt.show()


"""
Neutrino density and velocity perturbations
"""
fig, axs = plt.subplots(ncols = 2, layout = "constrained", figsize = (16, 6))

for i, k in enumerate(k_vals):
    axs[0].plot(data[k]["x"], 4*data[k]["Nu_0"], colors_nu[i], linewidth = 1)
    axs[1].plot(data[k]["x"], -3*data[k]["Nu_1"], colors_nu[i], label = r"$k=\,$" + f"{k}", linewidth = 1)
    axs[0].axvline(x_enter[i], color = colors_nu[i], linestyle = "--", linewidth = 1)
    axs[1].axvline(x_enter[i], color = colors_nu[i], linestyle = "--", linewidth = 1)
for i in range(2): #TODO: maybe no point in including here?
    axs[i].axvline(x_rm, color = "gold", label = r"$(r,m)$-eq." if i == 1 else None, linewidth = 1, path_effects = [pe.Stroke(linewidth = 2, foreground = "black"), pe.Normal()])
    axs[i].axvline(x_mLambda, color = "palevioletred", label = r"$(m,\Lambda)$-eq." if i == 1 else None, linewidth = 1, path_effects = [pe.Stroke(linewidth = 2, foreground = "black"), pe.Normal()])

axs[1].legend(loc = "lower left", framealpha = 1)
axs[0].set_xlim(x_min, x_max)
axs[1].set_xlim(x_min, x_max)
axs[0].set_title(r"Density $\delta_\nu=4\mathcal{N}_0$")
axs[1].set_title(r"Velocity $v_\nu=-3\mathcal{N}_1$")
fig.supxlabel(r"$x=\log(a)$")
fig.savefig("figs/neutrinos.pdf")
plt.show()


"""
Density and velocity perturbations
"""
fig, axs = plt.subplots(nrows = 3, ncols = 2, layout = "constrained", figsize = (16, 16))

lines = []
labels = []
for i, k in enumerate(k_vals):
    """
    CDM and baryons
    """
    l1 = axs[0][0].plot(data[k]["x"], data[k]["delta_b"], colors_b[i], linewidth = 1)[0]
    l2 = axs[0][0].plot(data[k]["x"], data[k]["delta_CDM"], colors_CDM[i], linewidth = 1)[0]
    lines.append((l1, l2))
    labels.append(r"$k=\,$" + f"{k}")
    axs[0][1].plot(data[k]["x"], data[k]["v_b"], colors_b[i], linewidth = 1)
    axs[0][1].plot(data[k]["x"], data[k]["v_CDM"], colors_CDM[i], linewidth = 1)
    for j in range(2):
        axs[0][j].axvline(x_enter[i], color = colors_b[i], linewidth = 1)
        axs[0][j].axvline(x_enter[i], color = colors_CDM[i], linestyle = "--", linewidth = 1)

    """
    Photons
    """
    axs[1][0].plot(data[k]["x"], 4*data[k]["Theta_0"], colors_gamma[i], label = r"$k=\,$" + f"{k}", linewidth = 1)
    axs[1][1].plot(data[k]["x"], -3*data[k]["Theta_1"], colors_gamma[i], linewidth = 1)
    axs[1][0].axvline(x_enter[i], color = colors_gamma[i], linestyle = "--", linewidth = 1)
    axs[1][1].axvline(x_enter[i], color = colors_gamma[i], linestyle = "--", linewidth = 1)

    """
    Neutrinos
    """
    axs[2][0].plot(data[k]["x"], 4*data[k]["Nu_0"], colors_nu[i], label = r"$k=\,$" + f"{k}", linewidth = 1)
    axs[2][1].plot(data[k]["x"], -3*data[k]["Nu_1"], colors_nu[i], linewidth = 1)
    axs[2][0].axvline(x_enter[i], color = colors_nu[i], linestyle = "--", linewidth = 1)
    axs[2][1].axvline(x_enter[i], color = colors_nu[i], linestyle = "--", linewidth = 1)

lines.append(axs[0][0].axvline(x_rm, color = "gold", linewidth = 1, path_effects = [pe.Stroke(linewidth = 2, foreground = "black"), pe.Normal()]))
lines.append(axs[0][0].axvline(x_mLambda, color = "palevioletred", linewidth = 1, path_effects = [pe.Stroke(linewidth = 2, foreground = "black"), pe.Normal()]))
labels.append(r"$(r,m)$-eq.")
labels.append(r"$(\Lambda,m)$-eq.")

titles = [[r"Matter densities $\delta_{\small\textrm{CDM}}$ (purple) and $\delta_b$ (orange)", r"Matter velocities $v_{\small\textrm{CDM}}$ (purple) and $v_b$ (orange)"],
          [r"Photon density $\delta_\gamma=4\Theta_0$", r"Photon velocity $v_\gamma=-3\Theta_1$"],
          [r"Neutrino density $\delta_\nu=4\mathcal{N}_0$", r"Neutrino velocity $v_\nu=-3\mathcal{N}_1$"]]
for i in range(3):
    for j in range(2):
        if not (i == 0 and j == 0):
            axs[i][j].axvline(x_rm, color = "gold", label = r"$(r,m)$-eq." if i == 0 else None, linewidth = 1, path_effects = [pe.Stroke(linewidth = 2, foreground = "black"), pe.Normal()])
            axs[i][j].axvline(x_mLambda, color = "palevioletred", label = r"$(m,\Lambda)$-eq." if i == 0 else None, linewidth = 1, path_effects = [pe.Stroke(linewidth = 2, foreground = "black"), pe.Normal()])
        axs[i][j].set_xlim(x_min, x_max)
        axs[i][j].set_title(titles[i][j])

axs[0][0].legend(lines, labels, handler_map = {tuple: HandlerTuple(ndivide = None)}, loc = "upper left", framealpha = 1)
axs[1][0].legend(loc = "lower left", framealpha = 1)
axs[2][0].legend(loc = "lower left", framealpha = 1)

axs[0][0].set_ylim(bottom = -8)
axs[0][0].set_yscale("asinh")
axs[0][1].set_yscale("asinh")
axs[0][0].set_yticks([1e4, 1e3, 1e2, 1e1, 2, 0, -2], labels = [r"$10^4$", r"$10^3$", "100", "10", "2", "0", "-2"])
axs[0][1].set_yticks([20, 10, 5, 2, 0.5, 0, -0.5, -2], labels = ["20", "10", "5", "2", "0.5", "0", "-0.5", "-2"])
fig.supxlabel(r"$x=\log(a)$")
fig.savefig("figs/density_and_velocity.pdf")
plt.show()


"""
Density perturbations
"""
fig, axs = plt.subplots(nrows = 2, ncols = 2, layout = "constrained", figsize = (16, 11))

quantities = [["delta_CDM", "delta_b"], ["Theta_0", "Nu_0"]]
factors = [[1, 1], [4, 4]]
colors = [[colors_CDM, colors_b], [colors_gamma, colors_nu]]
locs = [["upper left", "upper left"], ["lower left", "lower left"]]
titles = [[r"Cold dark matter $\delta_{\small\textrm{CDM}}$", r"Baryons $\delta_{b}$"], [r"Photons $\delta_\gamma=4\Theta_0$", r"Neutrinos $\delta_\nu=4\mathcal{N}_0$"]]
for i in range(2):
    for j in range(2):
        for k, k_val in enumerate(k_vals):
            axs[i][j].plot(data[k_val]["x"], factors[i][j]*data[k_val][quantities[i][j]], colors[i][j][k], label = r"$k=\,$" + f"{k_val}", linewidth = 1)
            axs[i][j].axvline(x_enter[k], color = colors[i][j][k], linestyle = "--", linewidth = 1)
        axs[i][j].axvline(x_rm, color = "gold", label = r"$(r,m)$-eq." if i == 0 and j == 0 else None, linewidth = 1, path_effects = [pe.Stroke(linewidth = 2, foreground = "black"), pe.Normal()])
        axs[i][j].axvline(x_mLambda, color = "palevioletred", label = r"$(m,\Lambda)$-eq." if i == 0 and j == 0 else None, linewidth = 1, path_effects = [pe.Stroke(linewidth = 2, foreground = "black"), pe.Normal()])
        axs[i][j].legend(loc = locs[i][j], framealpha = 1)
        axs[i][j].set_xlim(x_min, x_max)
        if i == 0:
            axs[i][j].set_yscale("asinh")
            axs[i][j].set_ylim(bottom = -8)
            axs[i][j].set_yticks([1e4, 1e3, 1e2, 1e1, 2, 0, -2], labels = [r"$10^4$", r"$10^3$", "100", "10", "2", "0", "-2"])
        axs[i][j].set_title(titles[i][j])

fig.supxlabel(r"$x=\log(a)$")
fig.savefig("figs/densities.pdf")
plt.show()


"""
Velocity perturbations
"""
fig, axs = plt.subplots(nrows = 2, ncols = 2, layout = "constrained", figsize = (16, 11))

quantities = [["v_CDM", "v_b"], ["Theta_1", "Nu_1"]]
factors = [[1, 1], [-3, -3]]
titles = [[r"Cold dark matter $v_{\small\textrm{CDM}}$", r"Baryons $v_{b}$"], [r"Photons $v_\gamma=-3\Theta_1$", r"Neutrinos $v_\nu=-3\mathcal{N}_1$"]]
for i in range(2):
    for j in range(2):
        for k, k_val in enumerate(k_vals):
            axs[i][j].plot(data[k_val]["x"], factors[i][j]*data[k_val][quantities[i][j]], colors[i][j][k], label = r"$k=\,$" + f"{k_val}", linewidth = 1)
            axs[i][j].axvline(x_enter[k], color = colors[i][j][k], linestyle = "--", linewidth = 1)
        axs[i][j].axvline(x_rm, color = "gold", label = r"$(r,m)$-eq." if i == 0 and j == 0 else None, linewidth = 1, path_effects = [pe.Stroke(linewidth = 2, foreground = "black"), pe.Normal()])
        axs[i][j].axvline(x_mLambda, color = "palevioletred", label = r"$(m,\Lambda)$-eq." if i == 0 and j == 0 else None, linewidth = 1, path_effects = [pe.Stroke(linewidth = 2, foreground = "black"), pe.Normal()])
        axs[i][j].legend(loc = locs[i][j], framealpha = 1)
        axs[i][j].set_xlim(x_min, x_max)
        if i == 0:
            axs[i][j].set_yscale("asinh")
            axs[i][j].set_yticks([20, 10, 5, 2, 0.5, 0, -0.5, -2], labels = ["20", "10", "5", "2", "0.5", "0", "-0.5", "-2"])
        axs[i][j].set_title(titles[i][j])

fig.supxlabel(r"$x=\log(a)$")
fig.savefig("figs/velocities.pdf")
plt.show()


"""
Photon and neutrino quadrupoles
"""
fig, axs = plt.subplots(ncols = 2, layout = "constrained", figsize = (16, 6))

lines = []
labels = []
for i, k in enumerate(k_vals):
    l1 = axs[0].plot(data[k]["x"], data[k]["Theta_2"], colors_gamma[i], linewidth = 1)[0]
    l2 = axs[1].plot(data[k]["x"], data[k]["Nu_2"], colors_nu[i], linewidth = 1)[0]
    lines.append((l1, l2))
    labels.append(r"$k=\,$" + f"{k}")
    axs[0].axvline(x_enter[i], color = colors_gamma[i], linestyle = "--", linewidth = 1)
    axs[1].axvline(x_enter[i], color = colors_nu[i], linestyle = "--", linewidth = 1)

axs[0].axvline(x_rm, color = "gold", linewidth = 1, path_effects = [pe.Stroke(linewidth = 2, foreground = "black"), pe.Normal()])
lines.append(axs[1].axvline(x_rm, color = "gold", linewidth = 1, path_effects = [pe.Stroke(linewidth = 2, foreground = "black"), pe.Normal()]))
labels.append(r"$(r,m)$-eq.")
axs[0].axvline(x_mLambda, color = "palevioletred", linewidth = 1, path_effects = [pe.Stroke(linewidth = 2, foreground = "black"), pe.Normal()])
lines.append(axs[1].axvline(x_mLambda, color = "palevioletred", linewidth = 1, path_effects = [pe.Stroke(linewidth = 2, foreground = "black"), pe.Normal()]))
labels.append(r"$(\Lambda,m)$-eq.")

axs[1].legend(lines, labels, handler_map = {tuple: HandlerTuple(ndivide = None)}, loc = "upper left", framealpha = 1)
axs[0].set_xlim(-12, x_max)
axs[1].set_xlim(x_min, x_max)
axs[0].set_title(r"Photon quadrupole $\Theta_2$")
axs[1].set_title(r"Neutrino quadrupole $\mathcal{N}_2$")
fig.supxlabel(r"$x=\log(a)$")
fig.savefig("figs/quadrupoles.pdf")
plt.show()


"""
Polarization multipoles
"""
fig, axs = plt.subplots(nrows = 3, layout = "constrained", figsize = (8, 16))

titles = [r"Monopole $\Theta^P_0$", r"Dipole $\Theta^P_1$", r"Quadrupole $\Theta^P_2$"]
for i, k in enumerate(k_vals):
    for j in range(3):
        axs[j].plot(data[k]["x"], data[k][f"Theta_P{j}"], colors_pol[i], label = r"$k=\,$" + f"{k}" if j == 1 else None, linewidth = 1)
        axs[j].axvline(x_enter[i], color = colors_pol[i], linestyle = "--", linewidth = 1)
        axs[j].set_title(titles[j])
        axs[j].set_xlim(-12, x_max)
for i in range(3):
    axs[i].axvline(x_rm, color = "gold", label = r"$(r,m)$-eq." if i == 2 else None, linewidth = 1, path_effects = [pe.Stroke(linewidth = 2, foreground = "black"), pe.Normal()])
    axs[i].axvline(x_mLambda, color = "palevioletred", label = r"$(m,\Lambda)$-eq." if i == 2 else None, linewidth = 1, path_effects = [pe.Stroke(linewidth = 2, foreground = "black"), pe.Normal()])

axs[1].legend(loc = "lower left", framealpha = 1)
axs[2].legend(loc = "lower left", framealpha = 1)
fig.supxlabel(r"$x=\log(a)$")
fig.savefig("figs/polarization.pdf")
plt.show()
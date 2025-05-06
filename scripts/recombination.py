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

# Reionization parameters
z_reion         = 8
Delta_z_reion   = 0.5
z_Hereion       = 3.5
Delta_z_Hereion = 0.5


"""
Define important time stamps
"""
x_Peebles = -7.50305
x_decoup_Saha = -7.15948
x_drag_Saha = -7.13972
x_recomb_Saha = -7.13113
x_decoup = -6.9666
x_drag = -6.94575
x_recomb = -6.94139
x_reion = np.log(1/(1 + z_reion))
x_Hereion = np.log(1/(1 + z_Hereion))


"""
Read data from files
"""
x_min = -12
x_max = 0
n_pts = (x_max - x_min)*10000 + 1

data = {}
keys = ["x", "Xe", "ne", "tau", "dtaudx", "ddtauddx", "g_tilde", "dgdx_tilde", "ddgddx_tilde", "s", "TCMB", "Tb", "tau_b", "dtau_bdx", "ddtau_bddx", "g_tilde_b", "dgdx_tilde_b", "ddgddx_tilde_b"]
for key in keys:
    data[key] = np.zeros(n_pts)

with open("results/recombination.txt", "r") as infile:
    for i, line in enumerate(infile):
        values = line.split()
        for j, key in enumerate(keys):
            data[key][i] = float(values[j])

data_Saha = {}
for key in keys[:3]:
    data_Saha[key] = np.zeros(n_pts)

with open("results/recombination_Saha.txt", "r") as infile:
    for i, line in enumerate(infile):
        values = line.split()
        for j, key in enumerate(keys[:3]):
            data_Saha[key][i] = float(values[j])


"""
Checking if g_tilde and g_tilde_b are properly normalized
"""
npts = 10000
integral = np.sum(data["g_tilde"])/npts
integral_b = np.sum(data["g_tilde_b"])/npts
rel_err = abs(integral - 1)
rel_err_b = abs(integral_b - 1)
print("    Photons:")
print(f"Integral of g_tilde: {integral}")
print(f"Relative error:      {rel_err*100:.3e}%")
print("    Baryons:")
print(f"Integral of g_tilde: {integral_b}")
print(f"Relative error:      {rel_err_b*100:.3e}%")


"""
Plot of free electron fraction X_e and free electron density n_e
"""
fig, axs = plt.subplots(ncols = 2, figsize = (16, 6))
fig.subplots_adjust(left = 0.04, right = 0.97, top = 0.9, wspace = 0.15)

Xe_ne = [data["Xe"], data["ne"], data_Saha["Xe"], data_Saha["ne"]]
for i in range(2):
    axs[i].plot(data["x"], Xe_ne[i+2], "k", linestyle = "--", label = r"Saha only", linewidth = 2)
    axs[i].plot(data["x"], Xe_ne[i], "k", label = r"Saha + Peebles", linewidth = 2)
    axs[i].axvline(x_Peebles, color = "slategrey", label = r"Saha$\,\to\,$Peebles")
    axs[i].axvline(x_recomb, color = "skyblue", label = r"Recombination")
    axs[i].axvline(x_recomb_Saha, color = "skyblue", linestyle = "--")
    axs[i].axvline(x_reion, color = "orange", label = r"Hydrogen reionization")
    axs[i].axvline(x_Hereion, color = "plum", label = r"Helium reionization")
    axs[i].set_xlim([x_min, x_max])

axs[1].legend(loc = "upper right", framealpha = 1)
axs[1].set_yscale("log")
axs[1].set_ylim(0.1*np.min(data["ne"]), 1.1*np.max(data["ne"]))
axs[0].set_title(r"Free electron fraction $X_e$")
axs[1].set_title(r"Electron number density $n_e$ [m$^{-3}$]")
fig.supxlabel(r"$x=\log(a)$")
fig.savefig("figs/Xe_and_ne.pdf")
plt.show()



"""
Plot of tau, -dtaudx and ddtauddx
"""
fig, axs = plt.subplots(ncols = 2, layout = "constrained", figsize = (16, 6))

quantities = [[data["tau"], -data["dtaudx"], data["ddtauddx"]], [data["tau_b"], -data["dtau_bdx"], data["ddtau_bddx"]]]
titles = ["Photon optical depth", "Baryon optical depth"]
for i in range(2): 
    axs[i].plot(data["x"], quantities[i][0], color = "cornflowerblue", label = r"$\tau(x)$" if i == 0 else None, linewidth = 2)
    axs[i].plot(data["x"], quantities[i][1], color = "palevioletred", label = r"$-\tau'(x)$" if i == 0 else None, linewidth = 2)
    axs[i].plot(data["x"], quantities[i][2], color = "yellowgreen", label = r"$\tau''(x)$" if i == 0 else None, linewidth = 2)
    axs[i].axvline(x_Peebles, color = "slategrey", label = r"Saha$\,\to\,$Peebles" if i == 1 else None)
    axs[i].axvline(x_decoup, color = "#fc95c7", label = "Photon decoupling" if i == 1 else None)
    axs[i].axvline(x_decoup_Saha, color = "#fc95c7", linestyle = "--")
    axs[i].axvline(x_drag, color = "gold", label = "Baryon decoupling" if i == 1 else None)
    axs[i].axvline(x_drag_Saha, color = "gold", linestyle = "--")
    axs[i].axvline(x_recomb, color = "skyblue", label = r"Recombination" if i == 1 else None)
    axs[i].axvline(x_recomb_Saha, color = "skyblue", linestyle = "--")
    axs[i].axvline(x_reion, color = "orange", label = "Hydrogen reionization" if i == 1 else None)
    axs[i].axvline(x_Hereion, color = "plum", label = "Helium reionization" if i == 1 else None)
    axs[i].legend(loc = "upper right", framealpha = 1)
    axs[i].set_title(titles[i])
    axs[i].set_xlim(x_min, x_max)
    axs[i].set_ylim(1e-6 if i == 0 else 1e-8, np.max(quantities[i][2]))
    axs[i].set_yscale("log") 

fig.supxlabel(r"$x=\log(a)$")
fig.supylabel(r"Optical depth and derivatives")
fig.savefig("figs/optical_depth.pdf")
plt.show()



# fig = plt.figure(layout = "constrained", figsize = (16, 12))
# subfigs = fig.subfigures(2, 1, hspace = 0.02, height_ratios = [1.5, 1])
# axs = list(subfigs[1].subplots(1, 2))
# axs.append(subfigs[0].subplots(1, 1))

# xlims = [[-9, -6], [-3.5, -0.5], [x_min, x_max]]
# ylims = [[1e-3, 1e5], [1e-7, 1], [1e-8, np.max(data["ddtau_bddx"])]]
# for i in range(3):
#     axs[i].plot(data["x"], data["tau"], color = "royalblue", label = r"$\tau(x)$", linewidth = 2, zorder = 1)
#     axs[i].plot(data["x"], -data["dtaudx"], color = "mediumvioletred", label = r"-$\tau'(x)$", linewidth = 2, zorder = 3)
#     axs[i].plot(data["x"], data["ddtauddx"], color = "olivedrab", label = r"$\tau''(x)$", linewidth = 2, zorder = 5)  
#     axs[i].axvline(x_Peebles, color = "slategrey", label = r"Saha$\,\to\,$Peebles", zorder = 6)
#     axs[i].axvline(x_decoup, color = "#fc95c7", label = "Photon decoupling", zorder = 7)
#     axs[i].axvline(x_decoup_Saha, color = "#fc95c7", linestyle = "--", zorder = 7)
#     axs[i].axvline(x_drag, color = "gold", label = "Baryon decoupling", zorder = 8)
#     axs[i].axvline(x_drag_Saha, color = "gold", linestyle = "--", zorder = 8)
#     axs[i].plot(data["x"], data["tau_b"], color = "cornflowerblue", label = r"$\tau_b(x)$", linewidth = 2, zorder = 0)
#     axs[i].plot(data["x"], -data["dtau_bdx"], color = "palevioletred", label = r"$-\tau_b'(x)$", linewidth = 2, zorder = 2)
#     axs[i].plot(data["x"], data["ddtau_bddx"], color = "yellowgreen", label = r"$\tau_b''(x)$", linewidth = 2, zorder = 4)
#     axs[i].axvline(x_recomb, color = "skyblue", label = r"Recombination")
#     axs[i].axvline(x_recomb_Saha, color = "skyblue", linestyle = "--")
#     axs[i].axvline(x_reion, color = "orange", label = "Hydrogen reionization")
#     axs[i].axvline(x_Hereion, color = "plum", label = "Helium reionization")
    
#     axs[i].set_xlim(xlims[i])
#     axs[i].set_ylim(ylims[i])
#     axs[i].set_yscale("log") 

# axs[2].legend(ncols = 2, framealpha = 1)
# fig.supxlabel(r"$x=\log(a)$")
# fig.supylabel(r"Optical depth and derivatives")
# fig.savefig("figs/optical_depth.pdf")
# plt.show()



"""
Plot of g_tilde, dgdx_tilde and ddgddx_tilde
"""
fig, axs = plt.subplots(nrows = 3, layout = "constrained", figsize = (8, 14))

quantities = [data["g_tilde"], data["dgdx_tilde"], data["ddgddx_tilde"]]
titles = [r"$\tilde{g}(x)$", r"$\tilde{g}'(x)$", r"$\tilde{g}''(x)$"]
for i in range(3):
    axs[i].plot(data["x"], quantities[i], "k", label = "Computed evolution", linewidth = 2, zorder = 1)
    axs[i].axvline(x_decoup, color = "#fc95c7", label = r"Photon decoupling") 
    axs[i].axvline(x_decoup_Saha, color = "#fc95c7", linestyle = "--")
    # axs[i].axvline(x_drag, color = "gold", label = r"Baryon decoupling") 
    # axs[i].axvline(x_drag_Saha, color = "gold", linestyle = "--") 
    axs[i].axvline(x_reion, color = "orange", label = "Hydrogen reionization")
    axs[i].axvline(x_Hereion, color = "plum", label = "Helium reionization")
    axs[i].set_title(titles[i])
    axs[i].set_xlim(-8 , x_max)
axs[0].legend(loc = "upper right", framealpha = 1)

fig.supxlabel(r"$x=\log(a)$")
fig.supylabel(r"Visibility function and derivatives")
fig.savefig("figs/visibility_function.pdf")
plt.show()




"""
Plot of sound horizon
"""
plt.figure(figsize = (8, 6))
plt.subplots_adjust(left = 0.16, right = 0.96, top = 0.93)
#TODO have important time stamps in this figure?
#TODO overplot close-ups of s(x), tau(x) and tau_b(x) around x_decoup and x_drag? 
plt.plot(data["x"], data["s"]/Mpc, "k", label = "Computed evolution", linewidth = 2)
# plt.axvline(x_Peebles, color = "slategrey", label = r"Saha$\,\to\,$Peebles") #TODO include this?
# plt.axvline(x_decoup, color = "#fc95c7", label = r"Photon decoupling") 
# plt.axvline(x_decoup_Saha, color = "#fc95c7", linestyle = "--") 
# plt.axvline(x_drag, color = "gold", label = r"Baryon decoupling") 
# plt.axvline(x_drag_Saha, color = "gold", linestyle = "--") 
plt.axvline(x_recomb, color = "skyblue", label = r"Recombination") 
plt.axvline(x_recomb_Saha, color = "skyblue", linestyle = "--") 
# plt.axvline(x_reion, color = "orange", label = "Hydrogen reionization") #TODO include this?
# plt.axvline(x_Hereion, color = "plum", label = "Helium reionization") #TODO include this?
plt.legend(loc = "lower right", framealpha = 1)
plt.xlim(x_min, x_max)
plt.yscale("log")
plt.yticks([1, 10, 100, 1000], ["1 Mpc", "10 Mpc", "100 Mpc", "1 Gpc"])
plt.xlabel(r"$x=\log(a)$")
plt.ylabel(r"Sound horizon $s(x)$")
plt.savefig("figs/sound_horizon.pdf")
plt.show()



"""
Plot of photon and baryon temperature versus redshift TODO decide if x or z
"""
# Temperatures today
print("\n")
print(f"T_CMB0 = {data['TCMB'][-1]:7.5f} K")
print(f"T_b0   = {data['Tb'][-1]*1e3:.5f} mK")

z = 1/np.exp(data["x"]) - 1

#TODO have important time stamps in this figure?
plt.figure(figsize = (8, 6))
plt.subplots_adjust(left = 0.12, right = 0.96, top = 0.93)
# plt.plot(z, data["TCMB"], color = "yellowgreen", label = r"$T_{\small\textrm{CMB}}$")
# plt.plot(z, data["Tb"], color = "#fca32d", label = r"$T_b$")
# plt.plot(data["x"], data["TCMB"], color = "yellowgreen", label = r"$T_{\small\textrm{CMB}}$")
# plt.plot(data["x"], data["Tb"], color = "#fca32d", label = r"$T_b$")
plt.plot(data["x"], data["TCMB"], "slategrey", label = r"$T_{\small\textrm{CMB}}$", linewidth = 2)
plt.plot(data["x"], data["Tb"], "k", label = r"$T_b$", linewidth = 2)
plt.plot(data["x"], np.exp(-data["x"])/np.exp(-data["x"][-1])*data["TCMB"][-1], color = "#9db92c", linestyle = "dashdot", label = r"$a^{-1}$")
plt.plot(data["x"], np.exp(-2*data["x"])/np.exp(-2*data["x"][-1])*data["Tb"][-1], color = "#a66fb5", linestyle = "dashdot", label = r"$a^{-2}$")
# plt.axvline(x_Peebles, color = "slategrey", label = r"Saha$\,\to\,$Peebles") #TODO include this?
plt.axvline(x_decoup, color = "#fc95c7", label = r"Photon decoupling") 
plt.axvline(x_decoup_Saha, color = "#fc95c7", linestyle = "--") 
plt.axvline(x_drag, color = "gold", label = r"Baryon decoupling") 
plt.axvline(x_drag_Saha, color = "gold", linestyle = "--") 
# plt.axvline(x_recomb, color = "skyblue", label = r"Recombination") 
# plt.axvline(x_recomb_Saha, color = "skyblue", linestyle = "--") 
# plt.axvline(x_reion, color = "orange", label = "Hydrogen reionization") #TODO include this?
# plt.axvline(x_Hereion, color = "plum", label = "Helium reionization") #TODO include this?
plt.legend(framealpha = 1)
# plt.xlim([np.max(z), z[-2]])
# plt.xscale("log")
plt.xlim([x_min, x_max])
plt.ylim([1e-2, 1e6])
plt.yscale("log")
# plt.xlabel(r"Redshift $z$")
plt.xlabel(r"$x=\log(a)$")
plt.ylabel(r"Temperature [K]")
plt.savefig("figs/photon_baryon_temp.pdf")
plt.show()



#TODO maybe plot helium and hydrogen density together with electron density?

#TODO compute tau without reionization as well? why split into three regions?
#TODO why should tau_b be computed without reionization? does it look wrong?

#TODO plot without helium for comparison?
#TODO compute predicted decoupling and recombination times from the Saha equation (without hydrogen?). overplot this with dashed but same colors
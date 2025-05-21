import numpy as np
import healpy as hp
import matplotlib.pyplot as plt
from matplotlib.legend_handler import HandlerTuple
from mpl_toolkits.axes_grid1.inset_locator import inset_axes, mark_inset
from matplotlib.colors import LinearSegmentedColormap
plt.rcParams.update({'text.usetex': True, 'font.size': 18, 'font.family': 'serif', 'font.serif': 'Computer Modern Sans Serif', 'font.weight': 100, 'mathtext.fontset': 'cm', 'xtick.labelsize': 16, 'ytick.labelsize': 16})

"""
Define important quantities
"""
h    = 0.67
k_eq = 0.0103782        # [1/Mpc]
Mpc  = 3.08567758e22    # [m]

"""
Read code output from files
"""
ell_max = 2500 
n_ells  = ell_max-1

k_min1  = 0.00001/Mpc
k_max1  = 0.3/Mpc
npts_k1 = int(np.log(k_max1) - np.log(k_min1))*10000 + 1

k_min2  = 0.000001
k_max2  = 5.0
npts_k2 = int(np.log(k_max2) - np.log(k_min2))*10000 + 1

r_min   = 1
r_max   = 500
npts_r  = int(r_max - r_min)*10000 + 1

npts    = [npts_k1, n_ells, n_ells, n_ells, n_ells, n_ells, n_ells, n_ells, npts_k2, 1001, 1001, npts_r]

data      = {}
filenames = ["transfer_functions_ell", "C_ells", "C_ells_1000", "C_ell_SW", "C_ell_ISW", "C_ell_Doppler", "C_ell_polarization", "C_ell_master", "P_k", "C_of_theta", "C_of_theta_1000", "xi"]
ells      = [15, 100, 350, 850, 1350, 1850]
keys      = [["k", "k*eta_0", "ThetaT", "ThetaE", "Nu", "Psi"], ["ell", "C_TT", "C_TE", "C_EE", "C_nu", "C_Psi", "C_lensed"], ["ell", "C_TT", "C_lensed"], ["ell", "C_TT"], ["ell", "C_TT"], ["ell", "C_TT"], ["ell", "C_TT"], ["ell", "C_TT"], ["k", "P_k"], ["theta", "C", "C_lensed"], ["theta", "C", "C_lensed"], ["r", "xi"]]
# keys      = [["k", "k*eta_0", "ThetaT", "ThetaE", "Nu", "Psi"], ["ell", "C_TT", "C_Psi", "C_lensed"], ["ell", "C_TT", "C_lensed"], ["ell", "C_TT"], ["ell", "C_TT"], ["ell", "C_TT"], ["ell", "C_TT"], ["ell", "C_TT"], ["k", "P_k"], ["theta", "C", "C_lensed"], ["theta", "C", "C_lensed"], ["r", "xi"]]
# keys      = [["k", "k*eta_0", "ThetaT", "ThetaE", "Nu", "Psi"], ["ell", "C_TT", "C_TE", "C_EE", "C_nu", "C_Psi"], ["ell", "C_TT"], ["ell", "C_TT"], ["ell", "C_TT"], ["ell", "C_TT"], ["k", "P_k"], ["theta", "C", "C_lensed"], ["r", "xi"]]
# keys      = [["k", "k*eta_0", "ThetaT", "ThetaE", "Nu"], ["ell", "C_TT", "C_TE", "C_EE", "C_nu", "C_Psi", "C_lensed"], ["ell", "C_TT"], ["ell", "C_TT"], ["ell", "C_TT"], ["ell", "C_TT"], ["k", "P_k"], ["theta", "C", "C_lensed"], ["r", "xi"]]
for i, filename in enumerate(filenames):
    data[filename] = {}
    if i == 0:
        for ell in ells:
            data[filename][ell] = {}
            for key in keys[i]:
                data[filename][ell][key] = np.zeros(npts[i])
            
            with open(f"results/{filename}{ell}.txt", "r") as infile:
                for j, line in enumerate(infile):
                    values = line.split()
                    for k, key in enumerate(keys[i]):
                        data[filename][ell][key][j] = float(values[k])
    else:
        for key in keys[i]:
            data[filename][key] = np.zeros(npts[i])

        with open(f"results/{filename}.txt", "r") as infile:
            for j, line in enumerate(infile):
                values = line.split()
                for k, key in enumerate(keys[i]):
                    data[filename][key][j] = float(values[k])

"""
CMB power spectrum data
"""
low_ell_TT  = {}
high_ell_TT = {}
high_ell_TE = {}
high_ell_EE = {}
theory      = {}
keys_low    = ["ell", "C_ell", "err_up", "err_down"]
keys_high   = ["ell", "C_ell", "err_down", "err_up", "bestfit"]
keys_theory = ["ell", "TT", "TE", "EE", "BB", "PP"]
for i in range(6):
    if i < 5:
        if i < 4:
            low_ell_TT[keys_low[i]] = []
        high_ell_TT[keys_high[i]]   = [] 
        high_ell_TE[keys_high[i]]   = [] 
        high_ell_EE[keys_high[i]]   = [] 
    theory[keys_theory[i]]          = []

with open(f"data/planck_cell_low.txt", "r") as infile:
    infile.readline()
    for line in infile:
        values = line.split()
        for k, key in enumerate(keys_low):
            low_ell_TT[key].append(float(values[k]))

spectra = ["TT-binned_R3.01", "TE-binned_R3.02", "EE-binned_R3.02"]
CMB_dicts = [high_ell_TT, high_ell_TE, high_ell_EE]
for spectrum, dict in zip(spectra, CMB_dicts):
    with open(f"data/COM_PowerSpect_CMB-{spectrum}.txt", "r") as infile:
        infile.readline()
        for line in infile:
            values = line.split()
            for k, key in enumerate(keys_high):
                dict[key].append(float(values[k]))

with open(f"data/COM_PowerSpect_CMB-base-plikHM-TTTEEE-lowl-lowE-lensing-minimum-theory_R3.01.txt", "r") as infile:
    infile.readline()
    for line in infile:
        values = line.split()
        for k, key in enumerate(keys_theory):
            theory[key].append(float(values[k]))

# Fix normalization of EE-spectrum
for key in keys_high[1:]:
    high_ell_EE[key] = 1e5 * 2*np.pi * np.array(high_ell_EE[key]) / (np.array(high_ell_EE["ell"]) * (np.array(high_ell_EE["ell"])+1))
theory["EE"] = 1e5 * 2*np.pi * np.array(theory["EE"]) / (np.array(theory["ell"]) * (np.array(theory["ell"])+1)) 

# Compute errors in ell-direction
high_ell_data = [high_ell_TT, high_ell_TE, high_ell_EE]
for high_ell in high_ell_data:
    high_ell["err_left"]  = []
    high_ell["err_right"] = []
    high_ell["err_left"].append((high_ell["ell"][1] - high_ell["ell"][0])/2)
    high_ell["err_right"].append((high_ell["ell"][1] - high_ell["ell"][0])/2)
    for i in range(1, len(high_ell["ell"])-1):
        high_ell["err_left"].append((high_ell["ell"][i] - high_ell["ell"][i-1])/2)
        high_ell["err_right"].append((high_ell["ell"][i+1] - high_ell["ell"][i])/2)
    high_ell["err_left"].append((high_ell["ell"][-1] - high_ell["ell"][-2])/2)
    high_ell["err_right"].append((high_ell["ell"][-1] - high_ell["ell"][-2])/2)

"""
Matter power spectrum data
"""
WMAP_ACT      = {}
SDSS_galaxies = {}
Lyman_alpha   = {}
keys          = ["k", "P_k", "err"]
for key in keys:
    WMAP_ACT[key]      = []
    SDSS_galaxies[key] = []
    Lyman_alpha[key]   = []

filenames = ["wmap_act", "reid_DR7", "lyalpha"]
P_k_dicts = [WMAP_ACT, SDSS_galaxies, Lyman_alpha]
for filename, dict in zip(filenames, P_k_dicts):
    with open(f"data/{filename}.txt", "r") as infile:
        infile.readline()
        for line in infile:
            values = line.split()
            for k, key in enumerate(keys):
                dict[key].append(float(values[k]))

# Compute the error in the WMAP/ACT data and Lyman-alpha forest data #TODO from SDSS BOSS survey 
WMAP_ACT["err"] = np.array(WMAP_ACT["err"]) - np.array(WMAP_ACT["P_k"])
Lyman_alpha["err"] = np.array(Lyman_alpha["err"]) - np.array(Lyman_alpha["P_k"])


"""
The transfer function Theta_ell(k) and the integrand |Theta_ell(k)|^2/k
"""
fig, axs = plt.subplots(ncols = 2, layout = "constrained", figsize = (16, 6))
axins = inset_axes(axs[1], width = 4.0, height = 3.0, loc = "upper right", bbox_to_anchor = (0.99, 0.99), bbox_transform = axs[1].transAxes)

x = np.array(data["transfer_functions_ell"][ells[0]]["k*eta_0"])
k = np.array(data["transfer_functions_ell"][ells[0]]["k"])*Mpc
colors_transfer = ["#fe797b", "#ffb750", "#ffea56", "#8fe968", "#36cedc", "#a587ca"]
for i, ell in enumerate(ells):
    ThetaT = np.array(data["transfer_functions_ell"][ell]["ThetaT"])
    integrand = abs(ThetaT)**2 / k
    axs[0].plot(x, np.sqrt(ell*(ell+1))*ThetaT, color = colors_transfer[i], label = r"$\ell=\,$" + f"{ell}", linewidth = 0.5) #TODO: correct normalization? divide by sqrt(2*pi) as well?
    axs[1].plot(x, ell*(ell+1)*integrand, color = colors_transfer[i], linewidth = 0.5) #TODO: correct normalization? divide by 2*pi as well?
    axins.plot(x, ell*(ell+1)*integrand, color = colors_transfer[i], linewidth = 0.5)

axs[0].legend(ncols = 2, framealpha = 1, columnspacing = 0.8)
axs[0].set_xlim(x[0], x[-1])
axs[1].set_xlim(x[0], x[-1])
axins.set_xlim(x[0], 2900)
axins.set_ylim(-0.3, 3.8)
axs[0].set_title(r"Transfer function $\sqrt{\ell(\ell+1)}\Theta_\ell(k)$")
axs[1].set_title(r"Integrand $\ell(\ell+1)|\Theta_\ell(k)|^2/k$ [Mpc]")
fig.supxlabel(r"$k\eta_0$")
fig.savefig("figs/transfer_function.pdf")
plt.show()


"""
The polarization transfer function Theta^E_ell(k) and the integrand |Theta^E_ell(k)|^2/k
"""
fig, axs = plt.subplots(ncols = 2, layout = "constrained", figsize = (16, 6))
axins = inset_axes(axs[1], width = 4.2, height = 3.0, loc = "upper right", bbox_to_anchor = (0.97, 0.95), bbox_transform = axs[1].transAxes)

x = np.array(data["transfer_functions_ell"][ells[0]]["k*eta_0"])
k = np.array(data["transfer_functions_ell"][ells[0]]["k"])*Mpc
colors_transfer = ["#fe797b", "#ffb750", "#ffea56", "#8fe968", "#36cedc", "#a587ca"]
for i, ell in enumerate(ells):
    ThetaE = np.array(data["transfer_functions_ell"][ell]["ThetaE"])
    integrand = abs(ThetaE)**2 / k
    axs[0].plot(x, np.sqrt((ell+2)*(ell+1)*ell*(ell-1))*ThetaE, color = colors_transfer[i], label = r"$\ell=\,$" + f"{ell}", linewidth = 0.5) #TODO: correct normalization?
    axs[1].plot(x, (ell+2)*(ell+1)*ell*(ell-1)*integrand, color = colors_transfer[i], label = r"$\ell=\,$" + f"{ell}", linewidth = 0.5) 
    axins.plot(x, (ell+2)*(ell+1)*ell*(ell-1)*integrand, color = colors_transfer[i], linewidth = 0.5)

axs[0].legend(ncols = 2, framealpha = 1, columnspacing = 0.8)
axs[0].set_xlim(x[0], x[-1])
axs[1].set_xlim(x[0], x[-1])
axins.set_xlim(x[0], 2900)
axins.set_ylim(-1e-8, 2.5e-7)
axs[0].set_title(r"Transfer function $\sqrt{(\ell+2)(\ell+1)\ell(\ell-1)}\Theta^E_\ell(k)$")
axs[1].set_title(r"Integrand $(\ell+2)(\ell+1)\ell(\ell-1)|\Theta^E_\ell(k)|^2/k$ [Mpc]")
fig.supxlabel(r"$k\eta_0$")
fig.savefig("figs/transfer_function_E.pdf")
plt.show()


"""
The CMB temperature power spectrum with data
"""
std = np.sqrt(2 / (2*np.array(data["C_ells"]["ell"]) + 1)) * np.array(data["C_ells"]["C_TT"]) #TODO correct?

fig, (ax) = plt.subplots(layout = "constrained", figsize = (16, 8))
axins = inset_axes(ax, width = 5.5, height = 3.2, loc = "upper left", bbox_to_anchor = (0.05, 0.99), bbox_transform = ax.transAxes)

ax.plot(data["C_ells"]["ell"], data["C_ells"]["C_TT"], "k", linewidth = 2, label = "Full spectrum", zorder = 2)
ax.plot(data["C_ell_master"]["ell"], data["C_ell_master"]["C_TT"], color = "#ffaed7", linewidth = 2, label = "Simplified spectrum", zorder = 1) 
ax.fill_between(data["C_ells"]["ell"], data["C_ells"]["C_TT"]+std, data["C_ells"]["C_TT"]-std, color = "k", alpha = 0.1, zorder = 0) #TODO: maybe remove
# ax.plot(theory["ell"], theory["TT"], color = "slategrey", linewidth = 2, label = "Planck base spectrum", zorder = 0) #TODO: maybe remove

axins.plot(data["C_ells"]["ell"], data["C_ells"]["C_TT"], "k", linewidth = 2, zorder = 2)
axins.plot(data["C_ell_master"]["ell"], data["C_ell_master"]["C_TT"], color = "#ffaed7", linewidth = 2, zorder = 1)
axins.fill_between(data["C_ells"]["ell"], data["C_ells"]["C_TT"]+std, data["C_ells"]["C_TT"]-std, color = "k", alpha = 0.1, zorder = 0) #TODO: maybe remove
# axins.plot(theory["ell"], theory["TT"], color = "slategrey", linewidth = 2, zorder = 0) #TODO: maybe remove

ax.errorbar(low_ell_TT["ell"], low_ell_TT["C_ell"], yerr = (low_ell_TT["err_up"], low_ell_TT["err_down"]), fmt = 'x', color = "#e26928", elinewidth = 1.5, capsize = 3, ecolor = "#f29b29", label = r"Low-$\ell$ Planck data")
axins.errorbar(low_ell_TT["ell"], low_ell_TT["C_ell"], yerr = (low_ell_TT["err_up"], low_ell_TT["err_down"]), fmt = 'x', color = "#e26928", elinewidth = 1.5, capsize = 3, ecolor = "#f29b29")
ax.errorbar(high_ell_TT["ell"], high_ell_TT["C_ell"], yerr = (high_ell_TT["err_up"], high_ell_TT["err_down"]), xerr = (high_ell_TT["err_left"], high_ell_TT["err_right"]), fmt = 'x', color = "#a64bb8", elinewidth = 1.5, capsize = 3, ecolor = "#d581e6", label = r"High-$\ell$ Planck data")
axins.errorbar(high_ell_TT["ell"], high_ell_TT["C_ell"], yerr = (high_ell_TT["err_up"], high_ell_TT["err_down"]), xerr = (high_ell_TT["err_left"], high_ell_TT["err_right"]), fmt = 'x', color = "#a64bb8", elinewidth = 1.5, capsize = 3, ecolor = "#d581e6")

ax.legend(loc = "upper right", framealpha = 1)
ax.set_xlim(2, ell_max)
axins.set_xlim(650, ell_max)
axins.set_ylim(-20, 2800)
ax.set_xscale("log")
axins.set_xscale("log")
ax.set_xticks([10, 100, 1000], ["10", "100", "1000"])
axins.set_xticks([700, 1000, 1500, 2000, 2500], ["700", "1000", "1500", "2000", ""])
fig.supxlabel(r"Multipole $\ell$")
fig.supylabel(r"$\ell(\ell+1)C_\ell^{\small\textrm{TT}}/2\pi$ [$\mu\textrm{K}^2$]")
plt.savefig("figs/CMB.pdf")
plt.show()


"""
The contributions to the CMB temperature power spectrum
"""
plt.figure(layout = "constrained", figsize = (9, 6))
plt.plot(data["C_ells"]["ell"], data["C_ells"]["C_TT"], "k", linewidth = 2, label = "Full spectrum")
plt.plot(data["C_ell_SW"]["ell"], data["C_ell_SW"]["C_TT"], color = "#f5b227", linewidth = 2, label = "SW")
plt.plot(data["C_ell_ISW"]["ell"], data["C_ell_ISW"]["C_TT"], color = "yellowgreen", linewidth = 2, label = "ISW")
plt.plot(data["C_ell_Doppler"]["ell"], data["C_ell_Doppler"]["C_TT"], color = "cornflowerblue", linewidth = 2, label = "Doppler")
plt.plot(data["C_ell_polarization"]["ell"], data["C_ell_polarization"]["C_TT"], color = "palevioletred", linewidth = 2, label = "Polarization")
plt.legend(loc = "upper left", bbox_to_anchor = [0.555, 0.39], framealpha = 1)
plt.xlim(2, ell_max)
plt.xscale("log")
plt.yscale("log")
plt.xticks([10, 100, 1000], ["10", "100", "1000"])
plt.xlabel(r"Multipole $\ell$")
plt.ylabel(r"$\ell(\ell+1)C_\ell^{\small\textrm{TT}}/2\pi$ [$\mu\textrm{K}^2$]")
plt.savefig("figs/CMB_contributions.pdf")
plt.show()


"""
The TE and EE polarization power spectra
"""
fig, axs = plt.subplots(ncols = 2, layout = "constrained", figsize = (16, 6))

lines = []
labels = []
axs[0].plot(data["C_ells"]["ell"], data["C_ells"]["C_TE"], "k", linewidth = 2, zorder = 1)
# axs[0].plot(theory["ell"], theory["TE"], color = "slategrey", linewidth = 2, zorder = 0) #TODO: maybe remove
l1 = axs[0].errorbar(high_ell_TE["ell"], high_ell_TE["C_ell"], yerr = (high_ell_TE["err_up"], high_ell_TE["err_down"]), xerr = (high_ell_TE["err_left"], high_ell_TE["err_right"]), fmt = 'x', color = "#3363bd", elinewidth = 1.5, capsize = 3, ecolor = "cornflowerblue")
lines.append(axs[1].plot(data["C_ells"]["ell"], data["C_ells"]["C_EE"], "k", linewidth = 2, zorder = 1)[0])
# lines.append(axs[1].plot(theory["ell"], theory["EE"], color = "slategrey", linewidth = 2, zorder = 0)[0]) #TODO: maybe remove
labels.append("Computed spectra")
# labels.append("Planck base spectra")
l2 = axs[1].errorbar(high_ell_EE["ell"], high_ell_EE["C_ell"], yerr = (high_ell_EE["err_up"], high_ell_EE["err_down"]), xerr = (high_ell_EE["err_left"], high_ell_EE["err_right"]), fmt = 'x', color = "mediumvioletred", elinewidth = 1.5, capsize = 3, ecolor = "palevioletred")
lines.append((l1, l2))
labels.append("Planck data")
axs[1].legend(lines, labels, handler_map = {tuple: HandlerTuple(ndivide = None)}, framealpha = 1)

axs[0].set_xlim(2, ell_max)
axs[1].set_xlim(2, ell_max)
# axs[1].set_ylim(-7, 98) #TODO maybe remove
axs[0].set_title(r"TE polarization spectrum")
axs[1].set_title(r"EE polarization spectrum")
axs[0].set_ylabel(r"$\ell(\ell+1)C_\ell^{\small\textrm{TE}}/2\pi$ [$\mu\textrm{K}^2$]")
axs[1].set_ylabel(r"$C_\ell^{\small\textrm{EE}}$ [$10^{-5}\mu\textrm{K}^2$]")
fig.supxlabel(r"Multipole $\ell$")
fig.savefig("figs/polarization_spectra.pdf")
plt.show()


"""
The neutrino power spectrum 
"""
plt.figure(layout = "constrained", figsize = (8, 6))
plt.plot(data["C_ells"]["ell"], data["C_ells"]["C_nu"], color = "skyblue", linewidth = 2)
plt.xlim(2, ell_max)
plt.xscale("log")
plt.xticks([10, 100, 1000], ["10", "100", "1000"])
plt.xlabel(r"Multipole $\ell$")
plt.ylabel(r"$\ell(\ell+1)C_\ell^{\nu}/2\pi$ [$\mu\textrm{K}^2$]")
plt.savefig("figs/neutrino_spectrum.pdf")
plt.show()


"""
The CMB lensing potential power spectrum
"""
plt.figure(layout = "constrained", figsize = (8, 6))
plt.plot(data["C_ells"]["ell"], data["C_ells"]["C_Psi"], "k", linewidth = 2, label = "Computed spectrum")
plt.plot(theory["ell"], np.array(theory["PP"])*1e7, color = "#ffaed7", linewidth = 2, linestyle = "-.", label = "Planck base spectrum")
plt.legend()
plt.xlim(2, ell_max)
plt.xscale("log")
plt.xticks([10, 100, 1000], ["10", "100", "1000"])
plt.xlabel(r"Multipole $\ell$")
plt.ylabel(r"$\ell^2(\ell+1)^2C_\ell^\Psi/2\pi$ [$10^{-7}$]") 
plt.savefig("figs/lensing_spectrum.pdf")
plt.show()


# """
# The normal and lensed angular correlation functions 
# """
# rel_diff = (np.array(data["C_of_theta"]["C"])-np.array(data["C_of_theta"]["C_lensed"]))/np.abs(np.array(data["C_of_theta"]["C"]))
# fig, axs = plt.subplots(ncols = 2, layout = "constrained", figsize = (16, 6))
# axins = inset_axes(axs[0], width = 4, height = 2.5, loc = "center left", bbox_to_anchor = (0.25, 0.5), bbox_transform = axs[0].transAxes)
# mark_inset(axs[0], axins, loc1 = 1, loc2 = 4, fc = "none", ec = "0.5")
# axs[0].plot(data["C_of_theta"]["theta"]*180/np.pi, data["C_of_theta"]["C"], "k", linewidth = 3, label = r"Unlensed")
# axs[0].plot(data["C_of_theta"]["theta"]*180/np.pi, data["C_of_theta"]["C_lensed"], color = "#ffaed7", linewidth = 1, label = r"Lensed")
# axins.plot(data["C_of_theta"]["theta"]*180/np.pi, data["C_of_theta"]["C"], "k", linewidth = 3)
# axins.plot(data["C_of_theta"]["theta"]*180/np.pi, data["C_of_theta"]["C_lensed"], color = "#ffaed7", linewidth = 1)
# axs[1].plot(data["C_of_theta"]["theta"]*180/np.pi, rel_diff, "k", linewidth = 0.75)
# axs[0].legend(framealpha = 1)
# axs[0].set_xlim(0, 180)
# axins.set_xlim(0, 5)
# axins.set_ylim(2.15e-10, np.max(np.array(data["C_of_theta"]["C"]))) 
# axs[1].set_xlim(0, 180)
# axs[0].set_title(r"Angular correlation functions $C(\theta)$")
# axs[1].set_title(r"Relative difference $(C(\theta)-C^{\small\textrm{lensed}}(\theta))/|C(\theta)|$")
# fig.supxlabel(r"Angular separation $\theta$ [deg]")
# plt.savefig("figs/angular_correlation.pdf")
# plt.show()


# """
# The lensed CMB temperature power spectrum with data
# """
# fig, (ax) = plt.subplots(layout = "constrained", figsize = (16, 8))
# axins = inset_axes(ax, width = 5.5, height = 3.2, loc = "upper left", bbox_to_anchor = (0.05, 0.99), bbox_transform = ax.transAxes)
# ax.plot(data["C_ells"]["ell"], data["C_ells"]["C_TT"], "k", linewidth = 3, label = "Unlensed spectrum")
# ax.plot(data["C_ells"]["ell"], data["C_ells"]["C_lensed"], color = "#ffaed7", linewidth = 1, label = "Lensed spectrum")
# axins.plot(data["C_ells"]["ell"], data["C_ells"]["C_TT"], "k", linewidth = 3)
# axins.plot(data["C_ells"]["ell"], data["C_ells"]["C_lensed"], color = "#ffaed7", linewidth = 1)
# ax.errorbar(low_ell_TT["ell"], low_ell_TT["C_ell"], yerr = (low_ell_TT["err_up"], low_ell_TT["err_down"]), fmt = 'x', color = "#e26928", elinewidth = 1.5, capsize = 3, ecolor = "#f29b29", label = r"Low-$\ell$ Planck data")
# axins.errorbar(low_ell_TT["ell"], low_ell_TT["C_ell"], yerr = (low_ell_TT["err_up"], low_ell_TT["err_down"]), fmt = 'x', color = "#e26928", elinewidth = 1.5, capsize = 3, ecolor = "#f29b29")
# ax.errorbar(high_ell_TT["ell"], high_ell_TT["C_ell"], yerr = (high_ell_TT["err_up"], high_ell_TT["err_down"]), fmt = 'x', color = "#a64bb8", elinewidth = 1.5, capsize = 3, ecolor = "#d581e6", label = r"High-$\ell$ Planck data")
# axins.errorbar(high_ell_TT["ell"], high_ell_TT["C_ell"], yerr = (high_ell_TT["err_up"], high_ell_TT["err_down"]), fmt = 'x', color = "#a64bb8", elinewidth = 1.5, capsize = 3, ecolor = "#d581e6")
# ax.legend(framealpha = 1)
# ax.set_xlim(2, ell_max)
# axins.set_xlim(650, ell_max)
# axins.set_ylim(0, 2800)
# ax.set_xscale("log")
# axins.set_xscale("log")
# ax.set_xticks([10, 100, 1000], ["10", "100", "1000"])
# axins.set_xticks([700, 1000, 1500, 2000, 2500], ["700", "1000", "1500", "2000", ""])
# fig.supxlabel(r"Multipole $\ell$")
# fig.supylabel(r"$\ell(\ell+1)C_\ell^{\small\textrm{TT}}/2\pi$ [$\mu\textrm{K}^2$]")
# plt.savefig("figs/CMB_lensed.pdf")
# plt.show()


"""
Lensed vs. unlensed angular correlation functions and power spectra
"""
rel_diff_corr = (np.array(data["C_of_theta"]["C"])-np.array(data["C_of_theta"]["C_lensed"]))/np.abs(np.array(data["C_of_theta"]["C"]))
rel_diff_corr_1000 = (np.array(data["C_of_theta_1000"]["C"])-np.array(data["C_of_theta_1000"]["C_lensed"]))/np.abs(np.array(data["C_of_theta_1000"]["C"]))
rel_diff_spec = (np.array(data["C_ells"]["C_TT"])-np.array(data["C_ells"]["C_lensed"]))/np.array(data["C_ells"]["C_TT"])
rel_diff_spec_1000 = (np.array(data["C_ells_1000"]["C_TT"])-np.array(data["C_ells_1000"]["C_lensed"]))/np.array(data["C_ells_1000"]["C_TT"])

fig, axs = plt.subplots(ncols = 2, layout = "constrained", figsize = (16, 6))
ax2_corr = axs[0].twinx()
ax2_spec = axs[1].twinx()
axins_corr = inset_axes(ax2_corr, width = 4, height = 2, loc = "center left", bbox_to_anchor = (0.15, 0.45), bbox_transform = ax2_corr.transAxes)
axins_spec = inset_axes(axs[1], width = 3.7, height = 1.9, loc = "center left", bbox_to_anchor = (0.19, 0.7), bbox_transform = axs[1].transAxes)

lines = []
labels = ["Unlensed", "Lensed", "Relative difference"]
lines.append(axs[0].plot(data["C_of_theta"]["theta"]*180/np.pi, data["C_of_theta"]["C"], "k", linewidth = 3, zorder = 0)[0])
l1 = axs[0].plot(data["C_of_theta"]["theta"]*180/np.pi, data["C_of_theta"]["C_lensed"], color = "#f6bad0", linewidth = 2, zorder = 2)[0]
l2 = axs[0].plot(data["C_of_theta_1000"]["theta"]*180/np.pi, data["C_of_theta_1000"]["C_lensed"], color = "#c90e53", linewidth = 2, zorder = 1)[0]
lines.append((l1, l2))
ax2_corr.plot(data["C_of_theta"]["theta"]*180/np.pi, rel_diff_corr, color = "#e1e38f", linewidth = 1, zorder = 2)
ax2_corr.plot(data["C_of_theta_1000"]["theta"]*180/np.pi, rel_diff_corr_1000, color = "#7a8018", linewidth = 1, zorder = 1)
axins_corr.plot(data["C_of_theta"]["theta"]*180/np.pi, rel_diff_corr, color = "#e1e38f", linewidth = 1, zorder = 2)
axins_corr.plot(data["C_of_theta_1000"]["theta"]*180/np.pi, rel_diff_corr_1000, color = "#7a8018", linewidth = 1, zorder = 1)

axs[1].plot(data["C_ells"]["ell"], data["C_ells"]["C_TT"], "k", linewidth = 3, zorder = 0)
axs[1].plot(data["C_ells"]["ell"], data["C_ells"]["C_lensed"], color = "#f6bad0", linewidth = 2, zorder = 2)
axs[1].plot(data["C_ells_1000"]["ell"], data["C_ells_1000"]["C_lensed"], color = "#c90e53", linewidth = 2, zorder = 1)
l3 = ax2_spec.plot(data["C_ells"]["ell"], rel_diff_spec, color = "#e1e38f", linewidth = 2, zorder = 2)[0]
l4 = ax2_spec.plot(data["C_ells_1000"]["ell"], rel_diff_spec_1000, color = "#7a8018", linewidth = 2, zorder = 1)[0]
lines.append((l3, l4))
axins_spec.plot(data["C_ells"]["ell"], data["C_ells"]["C_TT"], "k", linewidth = 3, zorder = 0)
axins_spec.plot(data["C_ells"]["ell"], data["C_ells"]["C_lensed"], color = "#f6bad0", linewidth = 2, zorder = 2)
axins_spec.plot(data["C_ells_1000"]["ell"], data["C_ells_1000"]["C_lensed"], color = "#c90e53", linewidth = 2, zorder = 1)

axs[0].legend(lines, labels, handler_map = {tuple: HandlerTuple(ndivide = None)}, loc = "upper left", framealpha = 1)
axs[0].set_xlim(0, 180)
axins_corr.set_xlim(110, 180)
axins_corr.set_ylim(-0.8e-5, 1e-5)
axins_corr.set_yticks([-0.7e-5, 0, 0.7e-5])
axs[1].set_xlim(2, ell_max)
axins_spec.set_xlim(1000, 2000)
axins_spec.set_ylim(100, 1400)
# axins_spec.set_yticks([-0.7e-5, 0, 0.7e-5]) #TODO
axs[0].set_xlabel(r"Angular separation $\theta$ [deg]")
axs[1].set_xlabel(r"Multipole $\ell$")
axs[0].set_ylabel(r"$C(\theta)$")
ax2_corr.set_ylabel(r"$(C(\theta)-C^{\small\textrm{lensed}}(\theta))/|C(\theta)|$")
axs[1].set_ylabel(r"$\ell(\ell+1)C_\ell^{\small\textrm{TT}}/2\pi$ [$\mu\textrm{K}^2$]")
ax2_spec.set_ylabel(r"$(C_\ell^{\small\textrm{TT}}-C_\ell^{\small\textrm{TT,lensed}})/C_\ell^{\small\textrm{TT}}$")
axs[0].set_title(r"Angular correlation functions")
axs[1].set_title(r"Power spectra")
plt.savefig("figs/lensed_vs_unlensed.pdf")
plt.show()


"""
The matter power spectrum
"""
fig, (ax) = plt.subplots(layout = "constrained", figsize = (16, 8))
axins = inset_axes(ax, width = 7, height = 4, loc = "lower left", bbox_to_anchor = (0.04, 0.07), bbox_transform = ax.transAxes)
# plt.figure(layout = "constrained", figsize = (8, 6))
ax.plot(data["P_k"]["k"], data["P_k"]["P_k"], "k", linewidth = 2, label = "Computed spectrum")
axins.plot(data["P_k"]["k"], data["P_k"]["P_k"], "k", linewidth = 2)
ax.axvline(k_eq/h, color = "gold", label = r"Equality scale $k_{\small\textrm{eq}}$") 
ax.errorbar(WMAP_ACT["k"], WMAP_ACT["P_k"], yerr = WMAP_ACT["err"], fmt = 'x', color = "#2199ca", elinewidth = 1.5, capsize = 3, ecolor = "#68c1e5", label = "CMB (WMAP+ACT)")
ax.errorbar(SDSS_galaxies["k"], SDSS_galaxies["P_k"], yerr = SDSS_galaxies["err"], fmt = 'x', color = "#d767ad", elinewidth = 1.5, capsize = 3, ecolor = "#ffaed7", label = "SDSS Galaxies (DR7 LRG)")
axins.errorbar(SDSS_galaxies["k"], SDSS_galaxies["P_k"], yerr = SDSS_galaxies["err"], fmt = 'x', color = "#d767ad", elinewidth = 1.5, capsize = 3, ecolor = "#ffaed7")
ax.errorbar(Lyman_alpha["k"], Lyman_alpha["P_k"], yerr = Lyman_alpha["err"], fmt = 'x', color = "olivedrab", elinewidth = 1.5, capsize = 3, ecolor = "yellowgreen", label = r"Ly$\alpha$ forest (BOSS)") #TODO correct that it is BOSS?
ax.legend(loc = "upper right", framealpha = 1)
ax.set_xlim(2e-3, 4) #TODO back to right=1?
ax.set_ylim(1, 5e4) #TODO back to bottom=70?
axins.set_xlim(5e-2, 3e-1)
axins.set_ylim(8e2, 2.5e4)
ax.set_xscale("log")
ax.set_yscale("log")
axins.set_xscale("log")
axins.set_yscale("log")
axins.set_yscale("log")
axins.set_xticks([5e-2, 6e-2, 7e-2, 8e-2, 9e-2, 1e-1, 2e-1, 3e-1], labels = ["", "0.06", "", "", "", "0.1", "0.2", ""])
fig.supxlabel(r"Wavenumber $k$ [$h/\textrm{Mpc}$]")
fig.supylabel(r"$P(k)$ [$(\textrm{Mpc}/h)^3$]")
fig.savefig("figs/P_k.pdf")
plt.show()


# #TODO: does it make sense to plot individual components?
# fig, axs = plt.subplots(ncols = 2, layout = "constrained", figsize = (16, 6))
# axs[0].plot(data["P_k"]["k"], data["P_k"]["P_k"], "k", linewidth = 2, label = "Computed spectrum")
# axs[0].axvline(k_eq/h, color = "k", linestyle = "--", label = r"Equality scale $k_{\small\textrm{eq}}$") 
# axs[0].errorbar(SDSS_galaxies["k"], SDSS_galaxies["P_k"], yerr = SDSS_galaxies["err"], fmt = 'x', color = "olivedrab", elinewidth = 1.5, capsize = 3, ecolor = "yellowgreen", label = "SDSS Galaxies (DR7 LRG)")
# axs[0].errorbar(WMAP_ACT["k"], WMAP_ACT["P_k"], yerr = WMAP_ACT["err"], fmt = 'x', color = "#2199ca", elinewidth = 1.5, capsize = 3, ecolor = "#68c1e5", label = "CMB (WMAP+ACT)")
# axs[1].plot(data["P_k"]["k"], data["P_k"]["P_k"], "k", linewidth = 2)
# axs[1].plot(data["P_k"]["k"], data["P_k"]["contrib_CDM"], color = "#eb6395", label = "CDM")
# axs[1].plot(data["P_k"]["k"], data["P_k"]["contrib_b"], color = "#f7a305", linestyle = "--", label = "Baryons")
# axs[1].plot(data["P_k"]["k"], data["P_k"]["contrib_gamma"], color = "#bec427", label = "Photons")
# axs[1].plot(data["P_k"]["k"], data["P_k"]["contrib_nu"], color = "#759cd0", label = "Neutrinos")
# axs[1].axvline(k_eq/h, color = "k", linestyle = "--")
# axs[0].legend(framealpha = 1)
# axs[1].legend(framealpha = 1)
# axs[0].set_xlim(2e-3, 1)
# axs[0].set_ylim(70, 5e4)
# axs[1].set_xlim(1e-4, 1)
# axs[1].set_ylim(5e-14, 5e5)
# axs[0].set_xscale("log")
# axs[1].set_xscale("log")
# axs[0].set_yscale("log")
# axs[1].set_yscale("log")
# axs[0].set_title("Total matter power spectrum")
# axs[1].set_title("Contributions from different components")
# fig.supxlabel(r"Wavenumber $k$ [$h/\textrm{Mpc}$]")
# fig.supylabel(r"$P(k)$ [$(\textrm{Mpc}/h)^3$]")
# plt.savefig("figs/P_k_components.pdf")
# plt.show()


"""
The correlation function
"""
r_drag = 149.10*h
plt.figure(layout = "constrained", figsize = (8, 6))
plt.plot(data["xi"]["r"], data["xi"]["xi"], "k", linewidth = 2, label = r"$r^2\xi(r)$")
plt.axvline(r_drag, color = "gold", label = r"$r_{\small\textrm{drag}}$") 
plt.legend(framealpha = 1)
# plt.xlim(r_min*h, r_max*h) #TODO
plt.xlim(r_min*h, 198) #TODO
plt.xlabel(r"Comoving separation $r$ [Mpc$/h$]")
plt.ylabel(r"Correlation function $r^2\xi(r)$ [(Mpc$/h)^2$]")
plt.savefig("figs/correlation_function.pdf")
plt.show()


"""
Generate a CMB map using HEALPix
"""
cmap = LinearSegmentedColormap.from_list(
    "cmap", 
    [(0.0,  "#02029e"), 
     (0.1,  "#0035fe"), 
     (0.2,  "#007bfe"), 
     (0.3,  "#00b7fe"), 
     (0.4,  "#51e1f1"), 
     (0.47, "#bee7e1"), 
     (0.53, "#fed98f"),
     (0.6,  "#fec139"), 
     (0.7,  "#fe9c00"), 
     (0.8,  "#fe6900"), 
     (0.9,  "#c62f00"), 
     (1.0,  "#600000"), 
    ],
    N = 256
)

C_ell = np.zeros(ell_max + 1)
C_ell[2:] = 2*np.pi / (data["C_ells"]["ell"]*(data["C_ells"]["ell"] + 1)) * data["C_ells"]["C_TT"] 
CMB_map = hp.synfast(C_ell, nside = 8192, lmax = ell_max, new = True)

rms = np.std(CMB_map)
hp.mollview(CMB_map, title = None, unit = r"$\mu$K", min = -3*rms, max = 3*rms, cmap = cmap, xsize = 1000)
plt.savefig("figs/CMB_map.pdf")
plt.show()
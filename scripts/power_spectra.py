import numpy as np
import matplotlib.pyplot as plt
from matplotlib.legend_handler import HandlerTuple
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
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
ell_max = 2000
n_ells  = ell_max-1

k_min1  = 0.00005/Mpc
k_max1  = 0.3/Mpc
npts_k1 = int(np.log(k_max1) - np.log(k_min1))*10000 + 1

k_min2  = 0.00005
k_max2  = 1.0
npts_k2 = int(np.log(k_max2) - np.log(k_min2))*10000 + 1

r_min   = 1
r_max   = 500
npts_r  = int(np.log(r_max) - np.log(r_min))*10000 + 1

npts    = [npts_k1, n_ells, n_ells, n_ells, n_ells, n_ells, npts_k2, npts_r]

data      = {}
filenames = ["transfer_functions_ell", "C_ells", "C_ell_SW", "C_ell_ISW", "C_ell_Doppler", "C_ell_polarization", "P_k", "xi"]
ells      = [15, 100, 350, 850, 1350, 1850]
keys      = [["k", "k*eta_0", "ThetaT", "ThetaE", "Nu", "Psi"], ["ell", "C_TT", "C_TE", "C_EE", "C_nu", "C_Psi"], ["ell", "C_TT"], ["ell", "C_TT"], ["ell", "C_TT"], ["ell", "C_TT"], ["k", "P_k", "contrib_CDM", "contrib_b", "contrib_gamma", "contrib_nu"], ["r", "xi"]]
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
keys_low    = ["ell", "C_ell", "err_up", "err_down"]
keys_high   = ["ell", "C_ell", "err_down", "err_up", "bestfit"]
for i in range(5):
    if i < 4:
        low_ell_TT[keys_low[i]] = []
    high_ell_TT[keys_high[i]]   = [] #TODO: use bestfit instead?
    high_ell_TE[keys_high[i]]   = [] #TODO: use bestfit instead?
    high_ell_EE[keys_high[i]]   = [] #TODO: use bestfit instead?

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

# Fix normalization of EE-spectrum TODO: change back
for key in keys_high[1:]:
    high_ell_EE[key] = 1e5 * 2*np.pi * np.array(high_ell_EE[key]) / (np.array(high_ell_EE["ell"]) * (np.array(high_ell_EE["ell"])+1))
#     high_ell_EE[key] = 1e5 * np.array(high_ell_EE[key])

"""
Matter power spectrum data
"""
SDSS_galaxies = {}
WMAP_ACT      = {}
keys          = ["k", "P_k", "err"]
for key in keys:
    SDSS_galaxies[key] = []
    WMAP_ACT[key]      = []

filenames = ["reid_DR7", "wmap_act"]
P_k_dicts = [SDSS_galaxies, WMAP_ACT]
for filename, dict in zip(filenames, P_k_dicts):
    with open(f"data/{filename}.txt", "r") as infile:
        infile.readline()
        for line in infile:
            values = line.split()
            for k, key in enumerate(keys):
                dict[key].append(float(values[k]))

# Compute the error in the WMAP/ACT data #TODO: correct? what is P_upper?
WMAP_ACT["err"] = np.array(WMAP_ACT["err"]) - np.array(WMAP_ACT["P_k"])


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
    #TODO: correct to use k*eta_0 as x-axis for both?

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
axins = inset_axes(axs[1], width = 4.0, height = 3.0, loc = "upper right", bbox_to_anchor = (0.99, 0.99), bbox_transform = axs[1].transAxes)

x = np.array(data["transfer_functions_ell"][ells[0]]["k*eta_0"])
k = np.array(data["transfer_functions_ell"][ells[0]]["k"])*Mpc
colors_transfer = ["#fe797b", "#ffb750", "#ffea56", "#8fe968", "#36cedc", "#a587ca"]
for i, ell in enumerate(ells):
    ThetaE = np.array(data["transfer_functions_ell"][ell]["ThetaE"])
    integrand = abs(ThetaE)**2 / k
    axs[0].plot(x, np.sqrt((ell+2)*(ell+1)*ell*(ell-1))*integrand, color = colors_transfer[i], label = r"$\ell=\,$" + f"{ell}", linewidth = 0.5) #TODO: correct normalization?
    axs[1].plot(x, (ell+2)*(ell+1)*ell*(ell-1)*integrand, color = colors_transfer[i], label = r"$\ell=\,$" + f"{ell}", linewidth = 0.5) 
    axins.plot(x, (ell+2)*(ell+1)*ell*(ell-1)*integrand, color = colors_transfer[i], linewidth = 0.5)

axs[0].legend(ncols = 2, framealpha = 1, columnspacing = 0.8)
axs[0].set_xlim(x[0], x[-1])
axs[1].set_xlim(x[0], x[-1])
axins.set_xlim(x[0], 2900)
axins.set_ylim(-0.3, 3.8)
axs[0].set_title(r"Transfer function $\sqrt{(\ell+2)(\ell+1)\ell(\ell-1)}\Theta^E_\ell(k)$")
axs[1].set_title(r"Integrand $(\ell+2)(\ell+1)\ell(\ell-1)|\Theta^E_\ell(k)|^2/k$ [Mpc]")
fig.supxlabel(r"$k\eta_0$")
fig.savefig("figs/transfer_function_E.pdf")
plt.show()


"""
The CMB temperature power spectrum with data
"""
#TODO: what about error in ell-direction?
#TODO: looks shifted in the ell-direction?
#TODO: theory results look slightly different from plot in milestone IV description?
fig, (ax) = plt.subplots(layout = "constrained", figsize = (16, 8))
axins = inset_axes(ax, width = 5.5, height = 3.2, loc = "upper left", bbox_to_anchor = (0.05, 0.99), bbox_transform = ax.transAxes)
ax.plot(data["C_ells"]["ell"], data["C_ells"]["C_TT"], "k", linewidth = 2, label = "Computed spectrum")
axins.plot(data["C_ells"]["ell"], data["C_ells"]["C_TT"], "k", linewidth = 2)
ax.errorbar(low_ell_TT["ell"], low_ell_TT["C_ell"], yerr = (low_ell_TT["err_up"], low_ell_TT["err_down"]), fmt = 'x', color = "#e26928", elinewidth = 1.5, capsize = 3, ecolor = "#f29b29", label = r"Low-$\ell$ Planck data")
axins.errorbar(low_ell_TT["ell"], low_ell_TT["C_ell"], yerr = (low_ell_TT["err_up"], low_ell_TT["err_down"]), fmt = 'x', color = "#e26928", elinewidth = 1.5, capsize = 3, ecolor = "#f29b29")
ax.errorbar(high_ell_TT["ell"], high_ell_TT["C_ell"], yerr = (high_ell_TT["err_up"], high_ell_TT["err_down"]), fmt = 'x', color = "#a64bb8", elinewidth = 1.5, capsize = 3, ecolor = "#d581e6", label = r"High-$\ell$ Planck data")
axins.errorbar(high_ell_TT["ell"], high_ell_TT["C_ell"], yerr = (high_ell_TT["err_up"], high_ell_TT["err_down"]), fmt = 'x', color = "#a64bb8", elinewidth = 1.5, capsize = 3, ecolor = "#d581e6")
ax.legend(framealpha = 1)
ax.set_xlim(2, ell_max)
axins.set_xlim(400, ell_max)
axins.set_ylim(0, 3000)
ax.set_xscale("log")
axins.set_xscale("log")
ax.set_xticks([10, 100, 1000], ["10", "100", "1000"])
axins.set_xticks([400, 500, 600, 700, 1000, 1500, 2000], ["", "500", "", "700", "1000", "1500", ""])
fig.supxlabel(r"Multipole $\ell$")
fig.supylabel(r"$\ell(\ell+1)C_\ell^{\small\textrm{TT}}/2\pi$ [$\mu\textrm{K}^2$]")
plt.savefig("figs/CMB.pdf")
plt.show()


"""
The contributions to the CMB temperature power spectrum
"""
plt.figure(layout = "constrained", figsize = (8, 6))
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
#TODO: both look shifted in the ell-direction?
axs[0].plot(data["C_ells"]["ell"], data["C_ells"]["C_TE"], "k", linewidth = 2)
l1 = axs[0].errorbar(high_ell_TE["ell"], high_ell_TE["C_ell"], yerr = (high_ell_TE["err_up"], high_ell_TE["err_down"]), fmt = 'x', color = "#3363bd", elinewidth = 1.5, capsize = 3, ecolor = "cornflowerblue")
# lines.append(axs[1].plot(data["C_ells"]["ell"], 1e-5*np.array(data["C_ells"]["C_EE"])*np.array(data["C_ells"]["ell"])*(np.array(data["C_ells"]["ell"])+1)/(2*np.pi), "k", linewidth = 2)[0]) #TODO: remove
lines.append(axs[1].plot(data["C_ells"]["ell"], data["C_ells"]["C_EE"], "k", linewidth = 2)[0])
labels.append("Computed spectra")
l2 = axs[1].errorbar(high_ell_EE["ell"], high_ell_EE["C_ell"], yerr = (high_ell_EE["err_up"], high_ell_EE["err_down"]), fmt = 'x', color = "mediumvioletred", elinewidth = 1.5, capsize = 3, ecolor = "palevioletred")
lines.append((l1, l2))
labels.append("Planck data")
axs[1].legend(lines, labels, handler_map = {tuple: HandlerTuple(ndivide = None)}, framealpha = 1)

axs[0].set_xlim(2, ell_max)
axs[1].set_xlim(2, ell_max)
axs[0].set_title(r"TE polarization spectrum")
axs[1].set_title(r"EE polarization spectrum")
axs[0].set_ylabel(r"$\ell(\ell+1)C_\ell^{\small\textrm{TE}}/2\pi$ [$\mu\textrm{K}^2$]")
axs[1].set_ylabel(r"$C_\ell^{\small\textrm{EE}}$ [$10^{-5}\mu\textrm{K}^2$]")
fig.supxlabel(r"Multipole $\ell$")
fig.savefig("figs/polarization_spectra.pdf")
plt.show()


"""
The matter power spectrum
"""
#TODO: find lyman alpha data?
#TODO: does it make sense to plot 2pi/r_drag?
r_drag = 149.10*h
plt.figure(layout = "constrained", figsize = (8, 6))
plt.plot(data["P_k"]["k"], data["P_k"]["P_k"], "k", linewidth = 2, label = "Computed spectrum")
plt.axvline(k_eq/h, color = "#ffaed7", label = r"Equality scale $k_{\small\textrm{eq}}$") 
plt.axvline(2*np.pi/r_drag, color = "gold", label = r"BAO scale $2\pi/r_{\small\textrm{drag}}$")
plt.errorbar(WMAP_ACT["k"], WMAP_ACT["P_k"], yerr = WMAP_ACT["err"], fmt = 'x', color = "#2199ca", elinewidth = 1.5, capsize = 3, ecolor = "#68c1e5", label = "CMB (WMAP+ACT)")
plt.errorbar(SDSS_galaxies["k"], SDSS_galaxies["P_k"], yerr = SDSS_galaxies["err"], fmt = 'x', color = "olivedrab", elinewidth = 1.5, capsize = 3, ecolor = "yellowgreen", label = "SDSS Galaxies (DR7 LRG)")
plt.legend(framealpha = 1)
plt.xlim(2e-3, 1)
plt.ylim(70, 5e4)
plt.xscale("log")
plt.yscale("log")
plt.xlabel(r"Wavenumber $k$ [$h/\textrm{Mpc}$]")
plt.ylabel(r"$P(k)$ [$(\textrm{Mpc}/h)^3$]")
plt.savefig("figs/P_k.pdf")
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
The neutrino power spectrum 
"""
#TODO: is this correct? maybe it makes sense?
plt.figure(layout = "constrained", figsize = (8, 6))
plt.plot(data["C_ells"]["ell"], data["C_ells"]["C_nu"], "k", linewidth = 2)
plt.xlim(2, ell_max)
plt.xscale("log")
plt.xticks([10, 100, 1000], ["10", "100", "1000"])
plt.xlabel(r"Multipole $\ell$")
plt.ylabel(r"$\ell(\ell+1)C_\ell^\nu/2\pi$ [$\mu\textrm{K}^2$]")
plt.savefig("figs/neutrino_spectrum.pdf")
plt.show()


"""
The CMB lensing potential power spectrum
"""
#TODO: fix
#TODO: find data? use public Planck theory spectrum?
plt.figure(layout = "constrained", figsize = (8, 6))
plt.plot(data["C_ells"]["ell"], data["C_ells"]["C_Psi"]*1e7, "k", linewidth = 2)
plt.xlim(2, ell_max)
# plt.xscale("log")
# plt.xticks([10, 100, 1000], ["10", "100", "1000"])
plt.xlabel(r"Multipole $\ell$")
plt.ylabel(r"$\ell^2(\ell+1)^2C_\ell^\Psi/2\pi$ [$10^{-7}$]") 
plt.savefig("figs/lensing_spectrum.pdf")
plt.show()


"""
The correlation function
"""
#TODO: fix
#TODO: find data?
#TODO: plot r_s too?
plt.figure(layout = "constrained", figsize = (8, 6))
plt.plot(data["xi"]["r"], data["xi"]["xi"], "k", linewidth = 2, label = r"Computed function")
plt.axvline(r_drag, color = "gold", label = r"$r_{\small\textrm{drag}}$") 
plt.legend(framealpha = 1)
plt.xlim(10, 210) #TODO: maybe change
plt.ylim(-100, 200) #TODO: maybe change
plt.xscale("log")
plt.xlabel(r"Comoving separation $r$ [$h^{-1}$Mpc]")
plt.ylabel(r"Correlation function $r^2\xi(r)$") #TODO: units?
plt.savefig("figs/correlation_function.pdf")
plt.show()
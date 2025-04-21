import numpy as np
import matplotlib.pyplot as plt
plt.rcParams.update({'text.usetex': True, 'font.size': 18, 'font.family': 'serif', 'font.serif': 'Computer Modern Sans Serif', 'font.weight': 100, 'mathtext.fontset': 'cm', 'xtick.labelsize': 16, 'ytick.labelsize': 16})

"""
Read data from files TODO: do this like Synne
"""
k_min   = 0.00005
k_max   = 0.3
n_pts_k = int((k_max - k_min)*10000 + 1)

ell_max = 2000
n_ells  = ell_max-1

n_pts   = [n_pts_k, n_ells]

data      = {}
filenames = ["P_k", "C_ells"]
keys      = [["k", "P_k"], ["ell", "C_TT"]]
for i, filename in enumerate(filenames):
    data[filename] = {}
    for key in keys[i]:
        data[filename][key] = np.zeros(n_pts[i])

    with open(f"results/{filename}.txt", "r") as infile:
        for j, line in enumerate(infile):
            values = line.split()
            for k, key in enumerate(keys[i]):
                data[filename][key][j] = float(values[k])

plt.plot(data["P_k"]["k"], data["P_k"]["P_k"])
plt.xscale("log")
plt.yscale("log")
plt.show()

plt.plot(data["C_ells"]["ell"], data["C_ells"]["C_TT"])
plt.xscale("log")
plt.yscale("log")
plt.show()
import numpy as np
import matplotlib.pyplot as plt
plt.rcParams.update({'text.usetex': True, 'font.size': 18, 'font.family': 'serif', 'font.serif': 'Computer Modern Sans Serif', 'font.weight': 100, 'mathtext.fontset': 'cm', 'xtick.labelsize': 16, 'ytick.labelsize': 16})

"""
Read data from files
"""
x_min = -18
x_max = 0
n_pts = (x_max - x_min)*10000 + 1

data = {}
k_vals = [0.001, 0.01, 0.1]
keys = ["x", "Phi", "Psi", "Pi", "delta_cdm", "delta_b", "v_cdm", "v_b", "Theta_0", "Theta_1", "Theta_2"]
for k in k_vals:
    data[k] = {}
    for key in keys:
        data[k][key] = np.zeros(n_pts)

    with open(f"results/perturbations_k{k}.txt", "r") as infile:
        for i, line in enumerate(infile):
            values = line.split()
            for j, key in enumerate(keys):
                data[k][key][i] = float(values[j])

colors = ["blue", "orange", "green"]

for i, k in enumerate(k_vals):
    plt.plot(data[k]["x"], data[k]["delta_cdm"], color = colors[i])
    plt.plot(data[k]["x"], np.abs(data[k]["delta_b"]), color = colors[i], linestyle = "--")
plt.ylim(bottom = 1e-1)
plt.yscale("log")
plt.show()

for i, k in enumerate(k_vals):
    plt.plot(data[k]["x"], data[k]["v_cdm"], color = colors[i])
    plt.plot(data[k]["x"], np.abs(data[k]["v_b"]), color = colors[i], linestyle = "--")
plt.ylim(bottom = 1e-6)
plt.yscale("log")
plt.show()

for i, k in enumerate(k_vals):
    plt.plot(data[k]["x"], data[k]["Theta_0"], color = colors[i])
plt.show()

for i, k in enumerate(k_vals):
    plt.plot(data[k]["x"], data[k]["Theta_1"], color = colors[i])
plt.show()

for i, k in enumerate(k_vals):
    plt.plot(data[k]["x"], data[k]["Theta_2"], color = colors[i])
plt.show()

for i, k in enumerate(k_vals):
    plt.plot(data[k]["x"], data[k]["Phi"], color = colors[i])
plt.show()

for i, k in enumerate(k_vals):
    plt.plot(data[k]["x"], data[k]["Phi"] + data[k]["Psi"], color = colors[i])
plt.show()
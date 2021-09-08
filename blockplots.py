import numpy as np
from matplotlib import pyplot as plt
from blocklib import *

filename = "timestep_data.dat"

zvals = []
ts_data = []

with open(filename,"r") as infile:
    for line in infile:
        t = line.strip().split(" ")
        if len(t) > 0:
            if t[0] == "z=":
                zvals.append(float(t[1]))
            else:
                ts_data.append([])
                for token in t:
                    ts_data[-1].append(float(token))

ts_data = np.array(ts_data)

# some specific z indices to plot
spec_z = [4, 7, 12, 17]

fig, axs = plt.subplots(nrows=4,ncols=3)
fig.suptitle("Error on Mean (relative to naive EoM) vs Number of Block Averages")

for i in range(len(spec_z)):
    axs[i][0].set_ylabel("z = {}".format(zvals[spec_z[i]]))
    for j in range(3):
        mu, c, sig, dsig = block_from_nparray(ts_data[(3*spec_z[i])+j],verbose=False)
        print("z={}, comp {}: Mean {}, Max Err {}".format(zvals[spec_z[i]],j,mu,np.amax(sig)))
        axs[i][j].errorbar(np.arange(len(sig)),sig/sig[0],yerr=dsig/sig[0])
    print("\n")

plt.show(block=False)
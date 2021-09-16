import numpy as np
from matplotlib import pyplot as plt
from blocking import *

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

# some specific z indices to plot
spec_z = [4, 7, 11, 17]

fig, axs = plt.subplots(nrows=4,ncols=3)
fig.suptitle("Error on Mean (relative to naive EoM) vs Number of Block Averages")

for i in range(len(spec_z)):
    axs[i][0].set_ylabel("z = {}".format(zvals[spec_z[i]]))
    for j in range(3):
        mu, c, sig, dsig = block(ts_data[(6*spec_z[i])+j][:512],verbose=False)
        if not c:
            continue
        fit_factor = np.sqrt( (1.+np.exp(-1./c[0][0]))/(1.-np.exp(-1./c[0][0])) )
        print("z={}, comp {}: Mean {}, Max Err {}".format(zvals[spec_z[i]],j,mu,np.amax(sig)))
        print("    Correlation time (in # of steps) from fit: {}".format(c[0][0]))
        print("    Asymptotic Err on Mean from fit: {}\n".format(fit_factor*sig[0]))
        axs[i][j].errorbar(np.arange(len(sig)),sig/sig[0],yerr=dsig/sig[0])
        kvals = np.linspace(0.0,len(sig),num=100)
        axs[i][j].plot(kvals, theory(kvals,*c[0]), "k:")
        axs[i][j].plot([0,kvals[-1]],[fit_factor]*2, "r:")
    print("\n")

plt.show(block=False)

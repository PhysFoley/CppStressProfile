import numpy as np
from matplotlib import pyplot as plt
import blocking

filename = "timestep_data.dat"
outfilename = "profile_eom.dat"

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

with open(outfilename,"w") as outfile:
    for i in range(len(zvals)):
        outfile.write("z= {}\n".format(zvals[i]))
        for j in range(6):
            mu, c, sig, dsig = blocking.block(ts_data[(6*i)+j],verbose=False)
            if not c:
                # if fitting fails, try again using data run with len = 2^n
                naive = sig[0]
                cutoff = int(np.power(2,np.floor(np.log2(len(ts_data[(6*i)+j])))))
                mu, c, sig, dsig = blocking.block(ts_data[(6*i)+j][:cutoff],verbose=False)
                if not c:
                    print("Falling back on naive Error on Mean")
                    outfile.write("{} ".format(naive))
                    continue
            fit_factor = np.sqrt( (1.+np.exp(-1./c[0][0]))/(1.-np.exp(-1./c[0][0])) )
            outfile.write("{} ".format(fit_factor*sig[0]))
        outfile.write("\n")

plt.show(block=False)

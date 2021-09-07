import numpy as np
from scipy.interpolate import CubicSpline
from matplotlib import pyplot as plt
from scipy.integrate import simps

# put all the input files here if you broke the trajectory up
filenames = ["stress_profile.dat"] # e.g. ["file1.dat","file2.dat","file3.dat"]

# same as above but for standard error files
errfiles = ["profile_std.dat"]

zvals = []
lat_profile = []
Szz_profile = []
n_files = len(filenames)

# read in the first file and populate the lists
with open(filenames[0],"r") as infile:
    for line in infile:
        t = line.strip().split(" ")
        if len(t) > 0:
            if t[0] == "z=":
                zvals.append(float(t[1]))
                lat_profile.append([])
                Szz_profile.append(0.0)
            else:
                t = [float(x) for x in t]
                lat_profile[-1].append( (0.5*(t[0]+t[4]) - t[8])/float(n_files) )
                Szz_profile[-1] += t[8]/float(n_files)

# now read in the rest of the data files and add their contributions
for fname in filenames[1:]:
    with open(fname,"r") as infile:
        z_index = -1
        kbn_index = 0
        for line in infile:
            t = line.strip().split(" ")
            if len(t) > 0:
                if t[0] == "z=":
                    z_index += 1
                    kbn_index = 0
                else:
                    t = [float(x) for x in t]
                    lat_profile[z_index][kbn_index] += (0.5*(t[0]+t[4]) - t[8])/float(n_files)
                    Szz_profile[z_index] += t[8]/float(n_files)
                    kbn_index += 1

lat_profile = np.array(lat_profile).transpose()
Szz_profile = np.array(Szz_profile)

# total lateral stress profile is the sum of the three contributions
# (kinetic, bonded, non-bonded)
tot_profile = lat_profile.sum(axis=0)

# now, read in error bar files
profile_err = [np.zeros(len(zvals)) for i in range(n_files)]
z_err = [np.zeros(len(zvals)) for i in range(n_files)]

for i in range(n_files):
    with open(errfiles[i],"r") as file:
        z_index = -1
        for line in file:
            t = line.strip().split(" ")
            if len(t) > 0:
                if t[0] == "z=":
                    z_index += 1
                else:
                    t = [float(x) for x in t]
                    profile_err[i][z_index] = np.sqrt( 0.25*(t[0]**2 + t[4]**2) + t[8]**2 )
                    z_err[i][z_index] = t[8]

tot_profile_err = []
tot_z_err = []
# this loop calculates the complete error bars based on all input files
for i in range(len(zvals)):
    p_sumsq = 0.0
    z_sumsq = 0.0
    for j in range(n_files):
        p_sumsq += profile_err[j][i]**2
        z_sumsq += z_err[j][i]**2
    tot_profile_err.append(np.sqrt(p_sumsq)/n_files)
    tot_z_err.append(np.sqrt(z_sumsq)/n_files)

tot_profile_err = np.array(tot_profile_err)
tot_z_err = np.array(tot_z_err)

half = int(np.floor(len(zvals)/2))

# Make some nice smooth cubic splines
p_cs = CubicSpline(zvals, tot_profile)
p_pe_cs = CubicSpline(zvals, tot_profile + tot_profile_err)
p_me_cs = CubicSpline(zvals, tot_profile - tot_profile_err)

z_cs = CubicSpline(zvals, Szz_profile)
z_pe_cs = CubicSpline(zvals, Szz_profile + tot_z_err)
z_me_cs = CubicSpline(zvals, Szz_profile - tot_z_err)

# integrate the stress profile for upper and lower leaflets using simpson's rule
lower_tension = simps(tot_profile[:half+1], zvals[:half+1])
upper_tension = simps(tot_profile[half:], zvals[half:])
full_tension = simps(tot_profile, zvals)

print("Upper leaflet tension: {}".format(upper_tension))
print("Lower leaflet tension: {}".format(lower_tension))
print("Total membrane tension: {}".format(full_tension))

plt.axhline(linewidth=1, color="black")
plt.axvline(linewidth=1, linestyle=":", color="black")

shift = zvals[half]

# z values for plotting the spline
zz = np.linspace(zvals[0], zvals[-1], 1000)

# plot the main splines
plt.plot(zz-shift, p_cs(zz) , zz-shift, z_cs(zz) )

#shade in the error regions
plt.fill_between(zz-shift, p_me_cs(zz), p_pe_cs(zz), alpha=0.5)
plt.fill_between(zz-shift, z_me_cs(zz), z_pe_cs(zz), alpha=0.5)

plt.xlabel(r"$z$ (from mid-plane) $[\sigma]$")
plt.ylabel(r"Lateral Stress $\partial\Sigma/\partial z$ $[\varepsilon/\sigma^3]$")

# plot the raw data
#plt.errorbar(zvals, tot_profile, yerr=tot_profile_err, linestyle="none")
#plt.errorbar(zvals, Szz_profile, yerr=tot_z_err, linestyle="none")

plt.show()

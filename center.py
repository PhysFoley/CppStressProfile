################################################################################
# centertraj.py
# Samuel Foley
#
# Center a VTF trajectory so that the bilayer midplane
# remains at the middle z-position of the box. This
# is a REQUIRED pre-processing step before using the
# stress calculation code. Note: Box must NOT fluctuate
# in z-direction.
#
# NOTE: This script simply places the overall centroid of
# bead positions at the middle of the box. It does not account
# for lipid number or shape asymmetry and generally the midplane
# of an asymmetric membrane will be off-center after being
# processed. Proper zeroing of the mid-plane in these cases
# will require modifying this script to take into account
# finer details for the situation of interest.
#
# The algorithm checks to see if there are more than a
# certain threshold number (max_num) of lipid beads near the
# edge of the simulation box (edge_size) and, if so, pre-shifts
# the membrane to avoid coordinate wrapping issues.
#
# One REQUIRED command line arg: input trajectory filename
# One optional argument: output trajectory filename
#
# Example usage:
#     $python centertraj.py traj.vtf centered_traj.vtf
#
################################################################################

import numpy as np
import sys

# these two parameters can be tweaked if the script
# is failing to properly center a trajectory
edge_size = 1.0
max_num = 5

ifname = "trajectory.vtf"
ofname = "centered_trajectory.vtf"

if len(sys.argv) < 2:
    print("Missing command line argument: VTF input filename")
    quit()
elif len(sys.argv) == 2:
    ifname = sys.argv[1]
else:
    ifname = sys.argv[1]
    ofname = sys.argv[2]

def load(filename, start=0, stop=-1, interval=1):
    step = []
    Lz = 0.0
    with open(filename,"r") as f:
        s = -1 # current step index
        for l in f:
            if(s > stop and stop != -1):
                break
            t = l.split()
            if(len(t) > 0 and t[0] != "bond" and t[0] != "unitcell" and t[0] != "atom" and t[0] != "timestep" and s >= start):
                step[-1].append(np.array([float(t[1]),float(t[2]),float(t[3])]))
            elif(len(t) > 0 and t[0] == "unitcell" and not Lz):
                Lz = float(t[3])
            elif(len(t) > 0 and t[0] == "timestep"): #we've reached the beginning of a new time step
                s += 1
                if(s%interval == 0 and s >= start and (s <= stop or stop == -1)):
                    step.append([])
    return step,Lz

def near_boundary(p_coords,box_z):
    num_high = 0
    num_low = 0
    
    for p in p_coords:
        folded = p[2]%box_z
        if folded > (box_z-edge_size):
            num_high += 1
        elif folded < edge_size:
            num_low += 1

    if (num_high > max_num) or (num_low > max_num):
        return True
    else:
        return False

step, box_z = load(ifname)

com = []
total_shift = []

for s in step:
    total_shift.append(0.0)
    shift = 0.0
    
    if near_boundary(s,box_z):
        shift = box_z/3.0
    total_shift[-1] += shift
    
    s = [p + np.array([0.0,0.0,shift]) for p in s]

    com = 0.0
    for p in s:
        com += (p[2]%box_z)/len(s)

    total_shift[-1] += (box_z/2.0) - com

# the required shift has now been calculated for all frames
# now we need to write out the modified trajectory
with open(ifname,"r") as orig, open(ofname,"w") as new:
    index = -1 #step index to be incrememented
    for line in orig:
        t = line.strip().split()
        if (len(t) == 0) or (t[0] == "unitcell") or (t[0] == "atom") or (t[0] == "bond"):
            new.write(line) #these lines go into the new file unmodified
        elif t[0] == "timestep":
            index += 1
            new.write(line)
        else:
            znew = float(t[3]) + total_shift[index]
            newline = t[0] + " " + t[1] + " " + t[2] + " " + str(znew) + "\n"
            new.write(newline)

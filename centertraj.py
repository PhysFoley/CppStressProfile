################################################################################
# centertraj.py
# Samuel Foley
#
# NOTE: This script is specificaly for 4-Bead Cooke model
# trajectories, and will not operate correctly for others
#
# Center a VTF trajectory so that the bilayer midplane
# remains at the middle z-position of the box. This
# is a REQUIRED pre-processing step before using the
# stress calculation code. Note: Box must NOT fluctuate
# in z-direction.
#
# Average tail bead positions are calculated separately
# for each monolayer, and then those results averaged
# to make sure the proper midplane is acquired for
# asymmetric membranes.
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

# Load coordinate data
def load(filename, start=0, stop=-1, interval=1):

    step = []
    coords = [0., 0., 0., 0.]
    Lz = None
    
    with open(filename,"r") as f:
        s = -1 # current step number

        for l in f:
            if(s > stop and stop != -1):
                break
            t = l.split()
            if(len(t) > 0 and t[0] != "bond" and t[0] != "unitcell" and t[0] != "atom" and t[0] != "timestep" and s >= start):
                if(int(t[0])%4 == 0): #head bead
                    coords[0] = np.array([float(t[1]),float(t[2]),float(t[3])])
                if(int(t[0])%4 == 1): #middle1 bead
                    coords[1] = np.array([float(t[1]),float(t[2]),float(t[3])])
                if(int(t[0])%4 == 2): #middle2 bead
                    coords[2] = np.array([float(t[1]),float(t[2]),float(t[3])])
                if(int(t[0])%4 == 3): #tail bead
                    coords[3] = np.array([float(t[1]),float(t[2]),float(t[3])])
                    step[-1].append(np.copy(coords))
            elif(len(t) > 0 and t[0] == "unitcell" and not Lz):
                Lz = float(t[3])
            elif(len(t) > 0 and t[0] == "timestep"): #we've reached the beginning of a new time step
                s += 1
                if(s%interval == 0 and s >= start and (s <= stop or stop == -1)):
                    step.append([])
    return step, Lz

# Check if a lipid is in the top leaflet
def intop(lipid):
    tth = np.array(lipid[0] - lipid[3])
    if np.dot(tth, np.array([0,0,1])) > 0:
        return True
    else:
        return False

step, box_z = load(ifname)

com = []
total_shift = []

#membrane shift algorithm
for s in step:
    total_shift.append(0.0)
    comtop = 0.0
    combot = 0.0
    com = 0.0
    num_high = 0
    num_low = 0
    shift = 0.0
    topcount = 0
    botcount = 0
    
    for lipid in s:
        folded = lipid[3][2]%box_z
        if intop(lipid):
            comtop += folded
            topcount += 1
        else:
            combot += folded
            botcount += 1
        if folded > box_z-edge_size:
            num_high += 1
        elif folded < edge_size:
            num_low += 1
    
    if num_high > max_num:
        shift = -box_z/3.0
    elif num_low > max_num:
        shift = box_z/3.0
    
    total_shift[-1] += shift
    
    if shift != 0.0:
        comtop = 0.0
        combot = 0.0
        for lipid in s:
            folded = (lipid[3][2]+shift)%box_z
            if intop(lipid):
                comtop += folded
            else:
                combot += folded
    
    comtop /= topcount
    combot /= botcount
    
    com = 0.5*(comtop + combot)
    
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



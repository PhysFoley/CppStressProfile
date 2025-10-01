# Membrane Stress Profile Calculator
Numerical software to calculate stress tensor components of laterally homogeneous slabs of membrane from MD simulations of the Cooke lipid model.

## Requirements

For the main code, the only requirement is a C++17 compiler; there are no dependencies beyond the standard libraries.

The Python pre- and post-analysis files centertraj.py, blockeom.py, and plotresults.py depend on NumPy, SciPy, and MatPlotLib.

## Usage

The functions defining the forces between beads of different types and their parameters are specified in a file called "params.hpp". This repository provides three such parameter headers:
* `default_params.hpp` : [Original Cooke lipid model](https://doi.org/10.1063/1.2135785)
* `flipfix_params.hpp` : [Flipfixed 4-bead Cooke lipids](https://doi.org/10.1021/acs.jctc.0c00862)
* `tapared_params.hpp` : [Tapered flipfixed 4-bead Cooke lipids](https://doi.org/10.1063/5.0189771)
In all cases, the conventions for bead type numbering follow those given in the simulation templates found [here](https://github.com/PhysFoley/pymbtools). Re-name the desired file to `params.hpp` before compiling, or write your own custom file tailored to your simulation. If `params.hpp` is not found, the code will compile with `default_params.hpp` and emit a warning.

To compile the main stress profile code, simply invoke `make` to generate the `stresscalc` executable.

The code requires trajectories to be in the VTF format written out by [ESPResSo MD](https://espressomd.org). The VTF trajectory file must be pre-processed such that the membrane is always located in the center of the box in the z-direction; this is done with
`python centertraj.py trajectory.vtf`
You can optionally give an output filename which by default is `centered_trajectory.vtf`.

The stresscalc routine is invoked with a statement like

`./stresscalc -f centered_trajectory.vtf -T 1.4 -d 10.0 -b 50 -e 499 -i 2 -n 25 -c 8`

* `-T 1.4` indicates the simulation was run with kT=1.4 ε
* `-d 10.0` thickness (z dimension) of the region at the center of the box in which to calculate the lateral stress profile
* `-b 50` begin with step 50 from the trajectory file
* `-e 499` end with step 499
* `-i 2` analyze every 2nd step (interval 2)
* `-n 25` calculate the stress at 25 evenly spaced points in the slab
* `-c 8` parallelize across 8 processes (when omitted, the code can automatically infer parallelization from SLURM environment variables)

After stress profile calculation completes (storing data by default in `timestep_data.dat` and `stress_profile.dat`), run `blockeom.py` to calculate error of the mean using blocking, and then `plotresults.py` to plot the stress profile with standard error shading. plotresults.py will also report the tension in the individual monolayer leaflets.

The output file stress_profile.dat contains line sequences like the following:

    z= 14.8
    -0.00381398 0 0 0 -0.00381398 0 0 0 -0.00381398
    0.00142565 0 0 -0.000170827 0.00167245 0 -0.000478575 0.000606222 0.0127527
    -0.00123832 0 0 -0.000141984 -0.0010852 0 0.000337088 -0.000500576 -0.00864861

The first line gives the z value to which the following three lines of data correspond. The first of the three data lines contains the kinetic contribution to the components of the stress tensor, in this case calculated from the ideal gas contribution. The second data line contains contributions to the components of the stress tensor due to bonded interactions. The third data line contains the contributions to the stress tensor components due to non-bonded interactions. The total stress tensor is obtained for a particular z-value by summing these lines component-wise. The component ordering is as follows: xx xy xz yx yy yz zx zy zz.

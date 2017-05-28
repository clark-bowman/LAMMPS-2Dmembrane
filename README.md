# LAMMPS-2Dmembrane

## 2d_membrane.cpp
C++ script which generates a LAMMPS script `2d_membrane.in` and data file `2d_membrane.dat`. Parameters are read in from `params.txt`.

## params.txt
Variable parameters for simulation. Each parameter is the only entry on a line. Currently, the only two parameters are top and bottom neutral angles for lipids (in degrees).

## 2d_membrane.in
LAMMPS script modeling a 2D bilipid membrane in a periodic simulation box filled with DPD fluid. `2d_membrane.dat` must be in the same directory to run. This script runs with the version of LAMMPS compiled Aug. 10, 2015.

_Outputs_:
* `curvature.out` - Average membrane curvature over time
* `2dmembrane.lammpstrj` - LAMMPS trajectory file for visualization. Loadable in, e.g., [VMD](http://www.ks.uiuc.edu/Research/vmd/). Does not include fluid. Commented out by default

Each output has two commented header lines, followed by a series of data in the form `timestep data` on each line. Timesteps are *not* the same as time units (convert with the timestep set in the script, default 0.03).

## 2d_membrane.dat
[Data file](http://lammps.sandia.gov/doc/read_data.html) which stores the initial positions of the lipids in a bilayer membrane. Includes positions, bonds, and angles.

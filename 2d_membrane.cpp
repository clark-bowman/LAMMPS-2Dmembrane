/*
LAMMPS Script Maker - 2D Bilipid Membrane, DPD Fluid (clark_bowman@brown.edu)
For LAMMPS version Aug. 10, 2015
*/

#include <fstream>
#include <stdlib.h>
#include <math.h>
#define PI 3.14159265359

int main()
{
    // Static parameters
    int tail_length = 6;
    int num_heads = 1;
    int num_tails = 2;
    double channel_length = 96.;
    double channel_width = 64.;
    double membrane_length = 64.;
    double head_spacing = 0.8;
    double tail_spacing = 0.5;
    double fluid_density = 2.0;
    double beadsize_head = 2.5;
    double pairwise_hh = 15.0;
    double pairwise_ww = 25.0;
    double pairwise_tt = 15.0;
    double pairwise_ht = 50.0;
    double pairwise_hw = 35.0;
    double pairwise_tw = 75.0;
    double angle_coeff = 50.0;
    double temp = 0.1;
    double dissipative_scaling = 1.0;
    double neutral_angle, neutral_angle_bot;
    int firsthead1, firsthead2, lasthead1, lasthead2, firsttail1, firsttail2, lasttail1, lasttail2;
    using namespace std;

    // Read dynamic parameters from text file
    {
        ifstream file_in ("params.txt");
        file_in >> neutral_angle;
        file_in >> neutral_angle_bot;
    }

    // Generate membrane initial condition
    {
        int n_horiz = int(ceil(membrane_length / (head_spacing * max(num_tails,num_heads))));
        ofstream fileOut ("2d_membrane.dat");
        fileOut.precision(5);
        fileOut << "# Membrane data file" << endl << endl;

        // Header with number of atoms, bonds, angles, dihedrals, and simulation information (dihedrals are not used here)
        fileOut << 2 * n_horiz * (num_heads + num_tails * tail_length) << " atoms" << endl;
        fileOut << 2 * n_horiz * (num_tails * tail_length + num_heads - 1) << " bonds" << endl;
        fileOut << 2 * n_horiz * (num_tails * (tail_length - 1) + 1) << " angles" << endl << endl;
        fileOut << "3 atom types" << endl << "2 bond types" << endl;
        fileOut << "3 angle types" << endl << endl;
        fileOut << 0. << " " << channel_length << " xlo xhi" << endl;
        fileOut << 0. << " " << channel_width << " ylo yhi" << endl;
        fileOut << "-0.1 0.1 zlo zhi" << endl << endl;

        // Masses for each particle type
        fileOut << "Masses" << endl << endl;
        fileOut << "1 " << pow(beadsize_head, 3) << endl << "2 1.0" << endl << "3 1.0" << endl << endl;

        // Loop to create atoms
        fileOut << "Atoms" << endl << endl;

        int mol_no = 1;
        int atom_no = 1;

        // Looping over lipids...
        for (int i = 0; i < n_horiz; i++)
        {
            // Up to the max number of heads/tails...
            for (int hctr = 0; hctr < max(num_heads,num_tails); hctr++)
            {
                // Create the head atom (if necessary)
                if (hctr < num_heads)
                {
                    fileOut << atom_no << " " << mol_no << " 1 " << 0.5 * (channel_length - membrane_length) + double(i * max(num_heads,num_tails) + hctr + 0.5*max(double(num_tails-num_heads),0.)) / double(n_horiz * max(num_heads,num_tails) - 1) * membrane_length << " " << 0.5 * channel_width - (double(tail_length) + 0.5) * tail_spacing << " 0" << endl;

                    // Certain lipid heads are marked for curvature calculation
                    if (i == 4)
                    {
                        firsthead2 = atom_no;
                        if(hctr == 0)
                            firsthead1 = atom_no;
                    }
                    if (i == n_horiz - 5)
                    {
                        lasthead2 = atom_no;
                        if(hctr == 0)
                            lasthead1 = atom_no;
                    }
                    atom_no++;
                }
                // If enough tails haven't been made yet...
                if (hctr < num_tails)
                {
                    // Initialize the tail location
                    double myx = 0.5 * (channel_length - membrane_length) + double(i * max(num_heads,num_tails) + hctr) / double(n_horiz * max(num_heads,num_tails) - 1) * membrane_length;
                    double myy = 0.5 * channel_width - (double(tail_length) + 0.5) * tail_spacing;

                    // Create each tail atom, moving along the length of the tail in the Y direction
                    for (int j = 0; j < tail_length; j++)
                    {
                        myy += tail_spacing;
                        fileOut << atom_no << " " << mol_no << " 2 " << myx << " " << myy << " 0" << endl;
                        atom_no++;
                    }

                    // Certain lipid tails are marked for curvature calculation
                    if (i == 4)
                    {
                        firsttail2 = atom_no - 1;
                        if(hctr == 0)
                            firsttail1 = atom_no - 1;
                    }
                    if (i == n_horiz - 5)
                    {
                        lasttail2 = atom_no - 1;
                        if(hctr == 0)
                            lasttail1 = atom_no - 1;
                    }
                }
            }
            // Update atom count for number of atoms in one lipid, increment molecule count by 1
            mol_no++;

            // Now do the same thing for the opposite layer of the membrane:
            for (int hctr = 0; hctr < max(num_heads,num_tails); hctr++)
            {
                if (hctr < num_heads)
                {
                    fileOut << atom_no << " " << mol_no << " 1 " << 0.5 * (channel_length - membrane_length) + double(i * max(num_heads,num_tails) + hctr + 0.5*max(double(num_tails-num_heads),0.)) / double(n_horiz * max(num_heads,num_tails) - 1) * membrane_length << " " << 0.5 * channel_width + (double(tail_length) + 0.5) * tail_spacing << " 0" << endl;
                    atom_no++;
                }
                if (hctr < num_tails)
                {
                    double myx = 0.5 * (channel_length - membrane_length) + double(i * max(num_heads,num_tails) + hctr) / double(n_horiz * max(num_heads,num_tails) - 1) * membrane_length;
                    double myy = 0.5 * channel_width + (double(tail_length) + 0.5) * tail_spacing;
                    for (int j = 0; j < tail_length; j++)
                    {
                        myy -= tail_spacing;
                        fileOut << atom_no << " " << mol_no << " 2 " << myx << " " << myy << " 0" << endl;
                        atom_no++;
                    }
                }
            }
            mol_no++;
        }

        // Loop to create bonds
        fileOut << endl << "Bonds" << endl << endl;
        int bondno = 1;

        // For each lipid...
        for (int i = 0; i < n_horiz * 2; i++)
        {
            for (int k = 0; k < max(num_heads,num_tails); k++)
            {
                int recenthead;
                if (k < num_heads)
                    recenthead = i * (tail_length * num_tails + num_heads) + k + 1 + min(k, num_tails) * tail_length;

                // The heads are bonded
                if (k < num_heads - 1)
                {
                    fileOut << bondno << " 2 " << recenthead << " " << i * (tail_length * num_tails + num_heads) + k + 2 + min(k + 1, num_tails) * tail_length << endl;
                    bondno++;
                }
                if (k < num_tails)
                {
                    // Each tail is bonded to a head
                    fileOut << bondno << " 1 " << recenthead << " " << i * (tail_length * num_tails + num_heads) + min(k+1,num_heads) + 1 + min(k, num_tails) * tail_length << endl;
                    bondno++;

                    // Atoms within a tail are bonded
                    for (int j = 1; j < tail_length; j++)
                    {
                        fileOut << bondno << " 1 " << i * (tail_length * num_tails + num_heads) + min(k+1,num_heads) + j + min(k, num_tails) * tail_length << " " << i * (tail_length * num_tails + num_heads) + min(k+1,num_heads) + 1 + j + min(k, num_tails) * tail_length << endl;
                        bondno++;
                    }
                }
            }
        }

        // Loop to create angles
        fileOut << endl << "Angles" << endl << endl;
        int angno = 1;

        // For each lipid...
        for (int i = 0; i < n_horiz * 2; i++)
        {
            // Add angle between first two tails
            fileOut << angno << " " << (i % 2) + 2 << " " << i * (tail_length * num_tails + num_heads) + 2 << " " << i * (tail_length * num_tails + num_heads) + 1 << " " << i * (tail_length * num_tails + num_heads) + tail_length + min(num_heads, 2) + 1 << endl;
            angno++;

            // Add angles along the length of each filament
            for (int k = 0; k < num_tails; k++)
            {
                int recenthead;
                if (k < num_heads)
                    recenthead = i * (tail_length * num_tails + num_heads) + 1 + k * (1 + tail_length);
                fileOut << angno << " 1 " << recenthead << " " << i * (tail_length * num_tails + num_heads) + min(k + 1,num_heads) + 1 + min(k, num_tails) * tail_length << " " << i * (tail_length * num_tails + num_heads) + min(k + 1,num_heads) + 2 + min(k, num_tails) * tail_length << endl;
                angno++;
                for (int j = 1; j < tail_length - 1; j++)
                {
                    fileOut << angno << " 1 " << i * (tail_length * num_tails + num_heads) + min(k + 1,num_heads) + j + min(k, num_tails) * tail_length << " " << i * (tail_length * num_tails + num_heads) + min(k + 1,num_heads) + 1 + j + min(k, num_tails) * tail_length << " " << i * (tail_length * num_tails + num_heads) + min(k + 1,num_heads) + 2 + j + min(k, num_tails) * tail_length << endl;
                    angno++;
                }
            }
        }
    }

    // Create LAMMPS script
    {
        ofstream fileOut ("2d_membrane.in");
        fileOut.precision(5);

        fileOut << "# LAMMPS Script - 2D Bilipid Membrane, DPD Fluid (clark_bowman@brown.edu)" << endl;
        fileOut << "# LAMMPS version Aug. 10, 2015" << endl << endl << endl << endl << endl << endl;

        fileOut << "# SEC: This section initializes the simulation." << endl << endl;
        fileOut << "# 	Define units, atom style, log path, and neighbor settings; read configuration data for the membrane." << endl;
        fileOut << "# 	Configuration data for the membrane in 2d_membrane.dat must be generated separately." << endl;
        fileOut << "# 	Final line communicates ghost data and is necessary for DPD parallelizing. On some versions of LAMMPS, use instead `communicate single vel yes'" << endl << endl;
        fileOut << "dimension 2" << endl << endl;
        fileOut << "units lj" << endl;
        fileOut << "atom_style molecular" << endl;
        fileOut << "log 2d_membrane.log" << endl;
        fileOut << "read_data 2d_membrane.dat" << endl;
        fileOut << "neighbor 0.3 bin" << endl;
        fileOut << "neigh_modify delay 3" << endl;
        fileOut << "comm_modify vel yes" << endl << endl;

        fileOut << "# 	Define bond, angle, pairwise interactions." << endl;
        fileOut << "# 	For pairwise interactions: type 1 is lipid head, 2 is lipid tail, 3 is water." << endl;
        fileOut << "#   914522 is a temperature seed; change for different randomness." << endl << endl;
        fileOut << "bond_style harmonic" << endl;
        fileOut << "bond_coeff 1 64.0 " << tail_spacing << endl;
        fileOut << "bond_coeff 2 64.0 " << head_spacing << endl;
        fileOut << "angle_style cosine/delta" << endl;
        fileOut << "angle_coeff 1 " << angle_coeff << " 180.0" << endl;
        fileOut << "angle_coeff 2 " << angle_coeff << " " << neutral_angle << endl;
        fileOut << "angle_coeff 3 " << angle_coeff << " " << neutral_angle_bot << endl << endl;
        fileOut << "pair_style dpd " << temp << " " << beadsize_head << " 914522" << endl;
        fileOut << "special_bonds lj 0.5 1.0 1.0" << endl;
        fileOut << "pair_coeff 1 1 " << pairwise_hh << " " << 4.5 * dissipative_scaling << " " << beadsize_head << endl;
        fileOut << "pair_coeff 1 3 " << pairwise_hw << " " << 4.5 * dissipative_scaling << " " << 0.5 * (beadsize_head + 1.) << endl;
        fileOut << "pair_coeff 2 3 " << pairwise_tw << " " << 20.0 * dissipative_scaling << " " << 1. << endl;
        fileOut << "pair_coeff 1 2 " << pairwise_ht << " " << 9.0 * dissipative_scaling << " " << 0.5 * (beadsize_head + 1.) << endl;
        fileOut << "pair_coeff 3 3 " << pairwise_ww << " " << 4.5 * dissipative_scaling << " 1.0" << endl;
        fileOut << "pair_coeff 2 2 " << pairwise_tt << " " << 4.5 * dissipative_scaling << " " << 1. << endl << endl;

        fileOut << "# SEC: This section initializes the geometry." << endl << endl;
        fileOut << "#   Define safe regions where fluid may be placed. This is region outside the membrane initialization block." << endl << endl;
        fileOut << "region safe block " << 0.5 * (channel_length - membrane_length - 1) << " " << 0.5 * (channel_length + membrane_length + 1) << " " << 0.5 * channel_width - (double(tail_length) + 1) * tail_spacing << " " << 0.5 * channel_width + (double(tail_length) + 1) * tail_spacing << " -0.1 0.1 side out" << endl;
        fileOut << "lattice sq "<< fluid_density << endl;
        fileOut << "create_atoms 3 region safe"<< endl << endl;

        fileOut << "# SEC: This section defines LAMMPS groups and computes that will be used in the simulation." << endl << endl;
        fileOut << "# 	Groups are named indicatively of their membership." << endl;
        fileOut << "# 	Computes include the positions of marked lipids, which are used to compute curvature." << endl << endl;
        fileOut << "group fluid type 3" << endl;
        fileOut << "group membrane type 1 2" << endl;
        fileOut << "group head1 id " << firsthead1;
        if (firsthead1 != firsthead2)
            fileOut << " " << firsthead2;
        fileOut << endl;
        fileOut << "group head2 id " << lasthead1;
        if (lasthead1 != lasthead2)
            fileOut << " " << lasthead2;
        fileOut << endl;
        fileOut << "group tail1 id " << firsttail1;
        if (firsttail1 != firsttail2)
            fileOut << " " << firsttail2;
        fileOut << endl;
        fileOut << "group tail2 id " << lasttail1;
        if (lasttail1 != lasttail2)
            fileOut << " " << lasttail2;
        fileOut << endl;
        fileOut << "compute mx membrane property/atom xu" << endl;
        fileOut << "compute my membrane property/atom yu" << endl;
        fileOut << "compute h1x head1 reduce ave c_mx" << endl;
        fileOut << "compute h1y head1 reduce ave c_my" << endl;
        fileOut << "compute h2x head2 reduce ave c_mx" << endl;
        fileOut << "compute h2y head2 reduce ave c_my" << endl;
        fileOut << "compute t1x tail1 reduce ave c_mx" << endl;
        fileOut << "compute t1y tail1 reduce ave c_my" << endl;
        fileOut << "compute t2x tail2 reduce ave c_mx" << endl;
        fileOut << "compute t2y tail2 reduce ave c_my" << endl;
        fileOut << "variable ang_pre equal 'abs(atan2(c_t1y-c_h1y,c_t1x-c_h1x)-atan2(c_t2y-c_h2y,c_t2x-c_h2x))'" << endl;
        fileOut << "variable ang equal 'v_ang_pre + 2 * (v_ang_pre > PI) * (PI - v_ang_pre)'" << endl;
        fileOut << "variable dist2 equal '(c_t1x-c_t2x)^2+(c_t1y-c_t2y)^2'" << endl;
        fileOut << "variable curv equal '2*sin(v_ang*0.5)/sqrt(v_dist2)'" << endl << endl;

        fileOut << "# SEC: This section initializes the particle velocities." << endl << endl;
        fileOut << "velocity all create " << temp << " 23473" << endl << endl;

        fileOut << "# SEC: This section defines fixes to impose forces in the simulation." << endl << endl;
        fileOut << "# 	NVE integration." << endl;
        fileOut << "# 	Fix 0 prevents particles from moving out of the X-Y plane." << endl;
        fileOut << "# 	Fix 4 records the curvature." << endl << endl;

        fileOut << "# 	Simulation timestep (LJ time units)." << endl;
        fileOut << "timestep 0.03" << endl;
        fileOut << "thermo 10000" << endl;
        fileOut << "fix 0 all enforce2d" << endl;
        fileOut << "fix 1 all nve/limit 0.05" << endl << endl;

        fileOut << "# 	Equilibration phase, followed by data recording." << endl;
        fileOut << "run 20000" << endl;
        fileOut << "fix 4 all ave/time 1 100 100 v_curv file curvature.out" << endl;
        fileOut << "#dump 1 all atom 500 2dmembrane.lammpstrj" << endl;
        fileOut << "run 30000" << endl << endl;
    }
    return 0;
}

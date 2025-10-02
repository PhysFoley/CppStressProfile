#pragma once

#include <iostream>
#include <sstream>
#include <cmath>
#include "forces.hpp"

//=========================================================
// Change these constants as appriopriate
//=========================================================

//make sure this matches the number of distinct bead types!
const int NUM_BEAD_TYPES = 2;

//make sure this is larger than longest-range non-bonded interaction!
const double NB_CUTOFF = 3.0;

//=========================================================
// Any user-defined globals can go here
//=========================================================

double bhh = 0.95;
double btt = 1.0;
double bht = 0.95;
double wc = 1.6;
double eps = 1.0;

//=========================================================
// Particle Pair Forces (see forces.hpp for built-in functions)
//=========================================================

double f00(double r, int i, int j)
{
    return lj(r, bhh, bhh*std::pow(2.0,1.0/6.0), eps);
}

double f11(double r, int i, int j)
{
    return ljcos2(r, btt, btt*std::pow(2.0,1.0/6.0), wc, eps);
}

double f01(double r, int i, int j)
{
    return lj(r, bht, bht*std::pow(2.0,1.0/6.0), eps);
}

//=========================================================
// Parse user-defined command-line arguments
//=========================================================

void parse_user_args(int argc, char** argv)
{
    return; // no extra args to parse for this simple model
}

//=========================================================
// Define Non-bonded Interaction Matrix
//=========================================================

void setup_interaction_matrix(double (*forces[NUM_BEAD_TYPES][NUM_BEAD_TYPES])(double,int,int))
{
    forces[0][0] = f00;
    forces[0][1] = forces[1][0] = f01;
    forces[1][1] = f11;
}

//=========================================================
// Bonded Interactions
//=========================================================

// these are used in the force calculation and must always be defined
double k_fene = 30.0;
double rmax_fene = 1.5;
double k_bend = 10.0;

// force between bonded particles with ids id1 and id2, and bead types type1 and type2
double bond_force(double r, double id1, double id2, double type1, double type2)
{
    // If particle IDs differ by 1, that's a FENE bond
    if(std::abs(id1-id2) == 1)
    {
        return fene(r, k_fene, rmax_fene);
    }
    else // otherwise, it should be harmonic (our "bending" force)
    {
        return harmonic(r, k_bend);
    }
}

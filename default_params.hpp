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
// Any new user-defined globals can go here
//=========================================================

double bhh = 0.95;
double btt = 1.0;
double bht = 0.95;
double wc = 1.6;
double eps = 1.0;

//=========================================================
// Particle Pair Forces
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
// Define Interaction Matrix
//=========================================================

void setup_interaction_matrix(double (*forces[NUM_BEAD_TYPES][NUM_BEAD_TYPES])(double,int,int))
{
    forces[0][0] = f00;
    forces[0][1] = forces[1][0] = f01;
    forces[1][1] = f11;
}

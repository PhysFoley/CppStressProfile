#pragma once

#include <sstream>
#include <cmath>
#include "forces.hpp"

//=========================================================
// Change these constants as appriopriate
//=========================================================

//make sure this matches the number of distinct bead types!
const int NUM_BEAD_TYPES = 8;

//make sure this is larger than longest-range non-bonded interaction!
const double NB_CUTOFF = 3.5;

//=========================================================
// Any new user-defined globals can go here
//=========================================================

double wc = 1.6;
double eps = 1.0;

// lipid taper angles (provided by user in cmd args below)
double alpha_a = 0.0;
double alpha_b = 0.0;

// matrix of relative diameters (populated in setup_interaction_matrix below)
double bb[NUM_BEAD_TYPES][NUM_BEAD_TYPES];

//=========================================================
// Particle Pair Forces
//=========================================================

// force involving head beads / flip-fixed beads
double fh(double r, int type1, int type2)
{
    return lj(r, bb[type1][type2], bb[type1][type2]*std::pow(2.0,1.0/6.0), eps);
}

// cohesive tail force
double ft(double r, int type1, int type2)
{
    return ljcos2(r, bb[type1][type2], bb[type1][type2]*std::pow(2.0,1.0/6.0), wc, eps);
}

//=========================================================
// Parse user-defined command-line arguments
//=========================================================

void parse_user_args(int argc, char** argv)
{
    std::stringstream ss;

    // IMPORTANT: start from i=1, because arg 0 is just the executable
    for(int i=1; i<argc-1; i+=2)
    {
        if(std::string(argv[i]).compare("-a1") == 0)
        {
            ss << argv[i+1];
            ss >> alpha_a;
            ss.clear();
        }
        else if(std::string(argv[i]).compare("-a2") == 0)
        {
            ss << argv[i+1];
            ss >> alpha_b;
            ss.clear();
        }
    }
}

//=========================================================
// Define Interaction Matrix
//=========================================================

void setup_interaction_matrix(double (*forces[NUM_BEAD_TYPES][NUM_BEAD_TYPES])(double,int,int))
{   
    alpha_a = alpha_a * (M_PI/180.0); // convert to radians
    alpha_b = alpha_b * (M_PI/180.0);

    double br[NUM_BEAD_TYPES];

    // type a bead radii
    br[3] = 0.5 + (1.5*std::sin(alpha_a));
    br[2] = 0.5 + (0.5*std::sin(alpha_a));
    br[1] = 0.5 - (0.5*std::sin(alpha_a));
    br[0] = 0.5 - (1.5*std::sin(alpha_a)) - 0.025;

    // type b bead radii
    br[7] = 0.5 + (1.5*std::sin(alpha_b));
    br[6] = 0.5 + (0.5*std::sin(alpha_b));
    br[5] = 0.5 - (0.5*std::sin(alpha_b));
    br[4] = 0.5 - (1.5*std::sin(alpha_b)) - 0.025;

    // b-param matrix (cross-bead "diameters")
    for(int i = 0; i < NUM_BEAD_TYPES; i++)
    {
        for(int j = 0; j < NUM_BEAD_TYPES; j++)
        {
            bb[i][j] = br[i] + br[j];
        }
    }
    
    // populate interaction matrix
    for(int i = 0; i < NUM_BEAD_TYPES; i++)
    {
        for(int j = 0; j <= i; j++) // note j <= i
        {
            // head beads
            if((i==0)||(i==4)||(j==0)||(j==4))
            {
                forces[i][j] = forces[j][i] = fh;
            }
            // flip-fix
            else if( ((j==1)||(j==2)) && ((i==5)||(i==6)) )
            {
                forces[i][j] = forces[j][i] = fh;
            }
            // all others
            else
            {
                forces[i][j] = forces[j][i] = ft;
            }
        }
    }
}


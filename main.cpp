/******************************************************************************
* main.cpp
* Samuel Foley
*
* For a complete description of how to use this program, consult the readme.
*
* Numerical routine for calculating Cooke membrane lateral stress profiles from
* simulation trajectories in VTF format. The central algorithm for calculating
* the stress profile is based on the Irving-Kirkwood formalism, and closely
* follows the presentation in Allen & Tildesley (2nd edition, pgs. 448 & 449).
*
* Warning: I'm not responsible if this program ruins your files, breaks
* your computer, and burns your house down. Use at your own risk.
*
*******************************************************************************/

#include <iostream>
#include <cstdlib>
#include <fstream>
#include <cmath>
#include <vector>
#include <queue>
#include <string>
#include <sstream>
#include <algorithm>
#include "mat3.hpp"

const double pi = 3.1415926535897;

//=========================================================
// IMPORTANT: Change these constants as appriopriate
//=========================================================

//make sure this matches the number of distinct bead types!
const int NUM_BEAD_TYPES = 7;

//make sure this is larger than longest-range non-bonded interaction!
const double NB_CUTOFF = 3.0;

//=========================================================
// Derivatives of Interaction Potentials
//=========================================================

double lj(double r, double b, double rc, double eps)
{
    if(r <= rc)
    {
        return 24.0*eps*( (std::pow(b,6.0)/std::pow(r,7.0)) - (2.0*std::pow(b,12.0)/std::pow(r,13.0)) );
    }
    else
    {
        return 0.0;
    }
}

double cos2(double r, double rc, double wc, double eps)
{
    if(r <= rc)
    {
        return 0.0;
    }
    else if(r > rc + wc)
    {
        return 0.0;
    }
    else
    {
        return (pi*eps/(2.0*wc))*std::sin(pi*(r-rc)/wc);
    }
}

double fene(double r, double k, double rinf)
{
    return k*r/(1.0 - std::pow(r/rinf,2.0) );
}

double bend(double r, double k)
{
    return k*(r-4.0);
}

double ljcos2(double r, double b, double eps, double rc, double wc)
{
    return lj(r, b, rc, eps) + cos2(r, rc, wc, eps);
}

double f00(double r)
{
    if(r <= 0.95*std::pow(2.0,1.0/6.0) )
    {
        return 24.0*( (std::pow(0.95,6.0)/std::pow(r,7.0)) - (2.0*std::pow(0.95,12.0)/std::pow(r,13.0)) );
    }
    else
    {
        return 0.0;
    }
}

double f01(double r)
{
    if(r <= 0.95*std::pow(2.0,1.0/6.0) )
    {
        return 24.0*( (std::pow(0.95,6.0)/std::pow(r,7.0)) - (2.0*std::pow(0.95,12.0)/std::pow(r,13.0)) );
    }
    else
    {
        return 0.0;
    }
}

double f11(double r)
{
    return lj(r, 1.0, std::pow(2.0,1.0/6.0), 1.0) + cos2(r, std::pow(2.0,1.0/6.0), 1.6, 1.0);
}

double f23(double r)
{
    if(r <= std::pow(2.0,1.0/6.0) )
        {
                return 24.0*( 1.0/std::pow(r,7.0) - (2.0/std::pow(r,13.0)) );
        }
        else
        {
                return 0.0;
        }
}

//=========================================================
// Algorithm Utility Functions
//=========================================================

double fold(double r, double box_l)
{
    if(r >= 0.0)
    {
        return std::fmod(r, box_l);
    }
    else
    {
        return box_l + std::fmod(r, box_l);
    }
}

double mind(double dr, double box_l)
{
    return dr - box_l*std::floor( (dr/box_l) + 0.5);
}

//Note: min_dr stores return value
void mind_3d(double* ri, double* rj, double* box_l, double* min_dr)
{
    for(int k=0; k<3; k++)
    {
        min_dr[k] = mind(ri[k]-rj[k], box_l[k]);
    }
}

double theta(double x)
{
    return (x >= 0.0) ? 1.0 : 0.0;
}

//=========================================================
// Trajectory Data Load Functions
//=========================================================

void split_str(std::string input, std::string delim, std::vector<std::string>* result)
{
    std::size_t pos = 0;
    std::size_t prev = 0;

    pos = input.find(delim);
    result->push_back(input.substr(prev,pos) );
    prev = pos;

    while(pos != std::string::npos)
    {
        pos = input.find(delim, prev+1);
        result->push_back(input.substr(prev+1,pos-prev) );
        prev = pos;
    }
}

void load_data(std::string filename, std::vector<std::vector<double*> >* step, std::vector<int>* type,
               std::vector<std::vector<int> >* bonds, std::vector<double*>* boxes, int start=0, int stop=-1, int interval=1)
{
    std::cout << "Loading trajectory from " << filename << std::endl;

    int update_inter = 100;
    if(stop != -1)
    {
        update_inter = std::ceil(float(stop - start)/10.0);
    }
    std::ifstream infile;
    infile.open(filename);
    std::string line;
    std::vector<std::string> t;

    double* init_box;

    std::stringstream ss;
    int inttemp = 0;
    int inttemp2 = 0;
    std::vector<std::string> svtemp;
    std::string stemp;
    int s = -1;
    while( !infile.eof() && (s <= stop || stop == -1) )
    {
        std::getline(infile, line);
        t.clear();
        split_str(line, " ", &t);

        if(t.size() > 1)
        {
            if(t[0].compare("atom") == 0)
            {
                ss.clear();
                ss << t[7];
                ss >> inttemp;
                type->push_back(inttemp);
                bonds->push_back(std::vector<int>() );
            }
            else if(t[0].compare("bond") == 0)
            {
                svtemp.clear();
                split_str(t[1], ":", &svtemp);
                ss.clear();
                ss << svtemp[0];
                ss >> inttemp;
                ss.clear();
                ss << svtemp[1];
                ss >> inttemp2;
                bonds->at(inttemp).push_back(inttemp2);
                bonds->at(inttemp2).push_back(inttemp);
            }
            else if(t[0].compare("timestep") == 0)
            {
                s++;
                if(s%interval == 0 && s >= start && (s <= stop || stop == -1) )
                {
                    if( ((s-start)/interval)%update_inter == 0)
                    {
                        std::cout << "Loaded " << (s-start)/interval << " steps..." << std::endl;
                    }
                    step->push_back(std::vector<double*>() );
                    ss = std::stringstream();
                }
            }
            else if(t[0].compare("unitcell") == 0)
            {
                double* box_l = new double[3];

                ss.clear();
                ss << t[1];
                ss >> box_l[0];
                ss << t[2];
                ss >> box_l[1];
                ss << t[3];
                ss >> box_l[2];

                if(s == -1)
                {
                    init_box = box_l;
                }
                else if(s%interval == 0 && s >= start && (s <= stop || stop == -1) )
                {
                    boxes->push_back(box_l);
                }
                else
                {
                    delete[] box_l;
                }
            }
            else if(s >= start)
            {
                double* coord = new double[3];

                ss.clear();
                ss << t[1];
                ss >> coord[0];
                ss.clear();
                ss << t[2];
                ss >> coord[1];
                ss.clear();
                ss << t[3];
                ss >> coord[2];

                if(s%interval == 0 && s >= start && (s <= stop || stop == -1) )
                {
                    step->back().push_back(coord);
                }
                else
                {
                    delete[] coord;
                }
            }
        }
    }
    infile.close();

    if(boxes->size() == 0)
    {
        boxes->push_back(init_box);
    }

    std::cout << "Finished Loading File.\n" << std::endl;
}

int main(int argc, char** argv)
{
    //=========================================================
    // Define Interaction Matrix
    //=========================================================

    double (*forces[NUM_BEAD_TYPES][NUM_BEAD_TYPES])(double);

    forces[0][0] = f00;
    forces[1][1] = f11;
    forces[2][2] = f11;
    forces[3][3] = f11;
    forces[4][4] = f11;
    forces[5][5] = f00;
    forces[6][6] = f00;

    forces[0][1] = forces[1][0] = f01;
    forces[0][2] = forces[2][0] = f01;
    forces[0][3] = forces[3][0] = f01;
    forces[0][4] = forces[4][0] = f01;
    forces[0][5] = forces[5][0] = f00;
    forces[0][6] = forces[6][0] = f00;

    forces[1][2] = forces[2][1] = f11;
    forces[1][3] = forces[3][1] = f11;
    forces[1][4] = forces[4][1] = f11;
    forces[1][5] = forces[5][1] = f01;
    forces[1][6] = forces[6][1] = f01;

    forces[2][3] = forces[3][2] = f23;
    forces[2][4] = forces[4][2] = f11;
    forces[2][5] = forces[5][2] = f01;
    forces[2][6] = forces[6][2] = f01;

    forces[3][4] = forces[4][3] = f11;
    forces[3][5] = forces[5][3] = f01;
    forces[3][6] = forces[6][3] = f01;
    forces[4][5] = forces[5][4] = f01;
    forces[4][6] = forces[6][4] = f01;

    forces[5][6] = forces[6][5] = f00;

    //=========================================================
    // User modifications should not be necessary
    // beyond this point
    //=========================================================

    //initialize command line args with default values
    std::stringstream ss;
    std::string ifname("centered_trajectory.vtf");
    std::string ofname("stress_profile.dat");
    std::string tfname("timestep_data.dat");
    int start = 0;
    int stop = -1;
    int interval = 1;
    int n_zvals = 51;
    double thickness = 10.0;
    double kT = 1.4;

    //parse command line args
    for(int i=1; i<argc-1; i+=2)
    {
        if(std::string(argv[i]).compare("-f") == 0)
        {
            ifname = argv[i+1];
        }
        else if(std::string(argv[i]).compare("-o") == 0)
        {
            ofname = argv[i+1];
        }
        else if(std::string(argv[i]).compare("-b") == 0)
        {
            ss << argv[i+1];
            ss >> start;
            ss.clear();
        }
        else if(std::string(argv[i]).compare("-e") == 0)
        {
            ss << argv[i+1];
            ss >> stop;
            ss.clear();
        }
        else if(std::string(argv[i]).compare("-i") == 0)
        {
            ss << argv[i+1];
            ss >> interval;
            ss.clear();
        }
        else if(std::string(argv[i]).compare("-n") == 0)
        {
            ss << argv[i+1];
            ss >> n_zvals;
            ss.clear();
        }
        else if(std::string(argv[i]).compare("-d") == 0)
        {
            ss << argv[i+1];
            ss >> thickness;
            ss.clear();
        }
        else if(std::string(argv[i]).compare("-t") == 0)
        {
            ss << argv[i+1];
            ss >> tfname;
            ss.clear();
        }
        else if(std::string(argv[i]).compare("-T") == 0)
        {
            ss << argv[i+1];
            ss >> kT;
            ss.clear();
        }
    }

    std::vector<std::vector<double*> > step;
    std::vector<int> type;
    std::vector<std::vector<int> > bonds;
    std::vector<double*> boxes;

    std::cout << std::endl;

    bool NVT = false; //default assumption, to be checked later
    load_data(ifname, &step, &type, &bonds, &boxes, start, stop, interval);

    std::cout << "Loaded:" << std::endl;
    std::cout << "       " << type.size() << " particles" << std::endl;
    std::cout << "       " << step.size() << " timesteps" << std::endl;
    std::cout << "       " << boxes.size() << " box geometries" << std::endl;
    std::cout << std::endl << std::endl;

    if(boxes.size() == 0 || type.size() == 0 || step.size() == 0)
    {
        std::cout << "Error reading trajectory, failed to load "
                  << "particles/timesteps/boxes" << std::endl;

        return EXIT_FAILURE;
    }

    if(boxes.size() != step.size())
    {
        if(boxes.size() == 1)
        {
            std::cout << "Only one box geometry found." << std::endl;
            std::cout << "Assuming constant volume {" << boxes[0][0] << ", "
                      << boxes[0][1] << ", " << boxes[0][2] << "}\n\n";

            NVT = true;
            for(int i = 1; i < step.size(); i++)
            {
                boxes.push_back(boxes[0]);
            }
        }
        else
        {
            std::cout << "Error: number of boxes and timesteps do not match";
            std::cout << std::endl;
            return EXIT_FAILURE;
        }
    }
    else if(boxes[0][2] != boxes[1][2])
    {
        std::cout << "Error: Box should not fluctuate in z-direction";
        std::cout << std::endl;
        return EXIT_FAILURE;
    }

    double mean_x = 0.0;
    double mean_y = 0.0;
    for(int i=0; i < boxes.size(); i++)
    {
      mean_x += boxes[i][0];
      mean_y += boxes[i][1];
    }
    mean_x = mean_x / float(boxes.size());
    mean_y = mean_y / float(boxes.size());
    std::cout << "Mean Lx: " << mean_x << std::endl;
    std::cout << "Mean Ly: " << mean_y << std::endl;

    double space = thickness/(n_zvals-1);
    double* zvals = new double[n_zvals];
    for(int i = 0; i < n_zvals; i++)
    {
        zvals[i]= (boxes[0][2]/2.0) - (thickness/2.0) + (i*space);
    }
    std::cout << std::endl << "Calculating stress at " << n_zvals;
    std::cout << " evenly spaced z values between ";
    std::cout << (boxes[0][2]/2.0) - (thickness/2.0);
    std::cout << " and " << (boxes[0][2]/2.0) + (thickness/2.0);
    std::cout << std::endl << std::endl;

    Mat3* Sk = new Mat3[n_zvals];
    Mat3* Sb = new Mat3[n_zvals];
    Mat3* Sn = new Mat3[n_zvals];

    Mat3* tmp_Sk = new Mat3[n_zvals];
    Mat3* tmp_Sb = new Mat3[n_zvals];
    Mat3* tmp_Sn = new Mat3[n_zvals];

    //array for storing step values to be block averaged
    double*** S_vals = new double**[n_zvals];
    for(int i = 0; i < n_zvals; i++)
    {
        S_vals[i] = new double*[9];
        for(int j = 0; j < 9; j++)
        {
            S_vals[i][j] = new double[step.size()];
        }
    }

    double r, phi_p;
    double *ri, *rj, *rb, *rij, *rib, *lo_r, *hi_r;
    double temp; // for storing temporary values
    int zbin_ind_i, zbin_ind_j, zbin_ind_b;

    rij = new double[3];
    rib = new double[3];

    double Lz = boxes[0][2]; //must be constant!

    int n_chkpts = 30; //length of progress bar

    std::cout << "Progress:" << std::endl;
    std::cout << "start|";
    for(int i = 0; i < n_chkpts; i++) { std::cout << " "; }
    std::cout << "|end" << std::endl;
    std::cout << "     |";
    std::cout.flush();

    std::queue<int> checkpoints;
    for(int i = 0; i < n_chkpts; i++)
    {
        checkpoints.push(i*step.size()/n_chkpts);
    }

    for(int s = 0; s < step.size(); s++)
    {
        //store value of stress tensor before new calculations
        for(int m = 0; m < n_zvals; m++)
        {
            tmp_Sk[m] = Sk[m];
            tmp_Sb[m] = Sb[m];
            tmp_Sn[m] = Sn[m];
        }
        //update the progress bar
        if(s == checkpoints.front())
        {
            std::cout << "*";
            std::cout.flush();
            checkpoints.pop();
        }
        double A = boxes[s][0]*boxes[s][1];
        for(int i = 0; i < type.size(); i++) //loop over all particles i
        {
            ri = step[s][i];
            zbin_ind_i = int((fold(ri[2],Lz)-zvals[0]+(space/2.0))/space);
            temp = kT/(A*space);
            for(int k = 0; k < 3; k++) //only diagonal components of kinetic tensor (ideal gas pressure)
            {
                //make sure the particle is within the domain we're treating
                if(zbin_ind_i < n_zvals && zbin_ind_i >= 0)
                {
                    Sk[zbin_ind_i].set(k,k,Sk[zbin_ind_i].get(k,k)-temp);
                }
            }
            for(int j = 0; j < i; j++) //loop over pairs i,j (without double-counting, hence j < i)
            {
                rj = step[s][j];
                mind_3d(ri, rj, boxes[s], rij); //rij stores return value
                r = std::sqrt(rij[0]*rij[0] + rij[1]*rij[1] + rij[2]*rij[2]);
                if(r > NB_CUTOFF)
                { //out of interaction range
                    continue;
                }
                else if(r == 0.0)
                {
                    std::cout << "Divide by zero in non-bonded interactions!" << std::endl;
                    return EXIT_FAILURE;
                }
                phi_p = forces[type[i]][type[j]](r);
                zbin_ind_j = int((fold(rj[2],Lz)-zvals[0]+(space/2.0))/space);
                if(zbin_ind_i == zbin_ind_j)
                {
                    //make sure the particle is within the domain we're treating
                    if(zbin_ind_i < n_zvals && zbin_ind_i >= 0)
                    {
                        //loop over tensor components
                        for(int a = 0; a < 3; a++)
                        {
                            for(int b = 0; b <= a; b++)
                            {
                                temp = rij[a]*rij[b]*phi_p*theta((fold(ri[2],Lz)-zvals[zbin_ind_i])/rij[2])*theta((zvals[zbin_ind_j]-fold(rj[2],Lz))/rij[2])/(A*space*r);
                                Sn[zbin_ind_i].set(a,b,Sn[zbin_ind_i].get(a,b)+temp);
                            }
                        }
                    }
                }
                else
                {
                    if(zbin_ind_i < zbin_ind_j)
                    {
                        lo_r = ri;
                        hi_r = rj;
                    }
                    else
                    {
                        lo_r = rj;
                        hi_r = ri;
                    }
                    int lo_ind = std::min(zbin_ind_i,zbin_ind_j);
                    int hi_ind = std::max(zbin_ind_i,zbin_ind_j);

                    //see pp. 449, Fig 14.1 Allen & Tildesley 2nd Ed.
                    double lo_factor = (zvals[std::max(lo_ind,0)]+(space/2.0)-fold(lo_r[2],Lz))/space;
                    double hi_factor = (fold(hi_r[2],Lz)-zvals[std::min(hi_ind,n_zvals-1)]+(space/2.0))/space;

                    //loop over z slabs between particles i and j
                    for(int m = lo_ind; m <= hi_ind; m++)
                    {
                        //only calculate if we're within calculation domain
                        if(m < n_zvals && m >= 0)
                        {
                            //loop over tensor components
                            for(int a = 0; a < 3; a++)
                            {
                                for(int b = 0; b <= a; b++)
                                {
                                    temp = rij[a]*rij[b]*phi_p*theta((fold(ri[2],Lz)-zvals[m])/rij[2])*theta((zvals[m]-fold(rj[2],Lz))/rij[2])/(A*std::abs(rij[2])*r);
                                    if(m == lo_ind) { temp = temp * lo_factor; }
                                    if(m == hi_ind) { temp = temp * hi_factor; }
                                    Sn[m].set(a,b,Sn[m].get(a,b)+temp);
                                }
                            }
                        }
                    }
                }
            } //end non-bonded interactions
            for(int bo = 0; bo < bonds[i].size(); bo++) //begin bonded interactions
            {
                rb = step[s][bonds[i][bo]];
                mind_3d(ri, rb, boxes[s], rib);
                r = std::sqrt(rib[0]*rib[0] + rib[1]*rib[1] +rib[2]*rib[2]);
                if(r == 0.0)
                {
                    std::cout << "Divide by zero in bonded interactions!" << std::endl;
                    return EXIT_FAILURE;
                }
                if(std::abs(bonds[i][bo]-i) == 1)
                {
                    phi_p = fene(r, 30.0, 1.5);
                }
                else
                {
                    phi_p = bend(r, 10.0);
                }
                zbin_ind_b = int((fold(rb[2],Lz)-zvals[0]+(space/2.0))/space);
                if(zbin_ind_i == zbin_ind_b)
                {
                    //make sure the particle is within the domain we're treating
                    if(zbin_ind_i < n_zvals && zbin_ind_i >= 0)
                    {
                        //loop over tensor components
                        for(int a = 0; a < 3; a++)
                        {
                            for(int b = 0; b <= a; b++)
                            {
                                temp = rib[a]*rib[b]*phi_p*theta((fold(ri[2],Lz)-zvals[zbin_ind_i])/rib[2])*theta((zvals[zbin_ind_b]-fold(rb[2],Lz))/rib[2])/(A*space*r);
                                Sb[zbin_ind_i].set(a,b,Sb[zbin_ind_i].get(a,b)+temp);
                            }
                        }
                    }
                }
                else
                {
                    if(zbin_ind_i < zbin_ind_b)
                    {
                        lo_r = ri;
                        hi_r = rb;
                    }
                    else
                    {
                        lo_r = rb;
                        hi_r = ri;
                    }
                    int lo_ind = std::min(zbin_ind_i,zbin_ind_b);
                    int hi_ind = std::max(zbin_ind_i,zbin_ind_b);

                    //see pp. 449, Fig 14.1 Allen & Tildesley 2nd Ed.
                    double lo_factor = (zvals[std::max(lo_ind,0)]+(space/2.0)-fold(lo_r[2],Lz))/space;
                    double hi_factor = (fold(hi_r[2],Lz)-zvals[std::min(hi_ind,n_zvals-1)]+(space/2.0))/space;

                    //loop over z slabs between particles i and j
                    for(int m = lo_ind; m <= hi_ind; m++)
                    {
                        //only calculate if we're within calculation domain
                        if(m < n_zvals && m >= 0)
                        {
                            //loop over tensor components
                            for(int a = 0; a < 3; a++)
                            {
                                for(int b = 0; b <= a; b++)
                                {
                                    temp = 0.5*rib[a]*rib[b]*phi_p*theta((fold(ri[2],Lz)-zvals[m])/rib[2])*theta((zvals[m]-fold(rb[2],Lz))/rib[2])/(A*std::abs(rib[2])*r);
                                    if(m == lo_ind) { temp = temp * lo_factor; }
                                    if(m == hi_ind) { temp = temp * hi_factor; }
                                    Sb[m].set(a,b,Sb[m].get(a,b)+temp);
                                }
                            }
                        }
                    }
                }
            }//end bonded interactions
        }//end particle loop

        //get this step's tensor values and store them for blocking later
        for(int m = 0; m < n_zvals; m++)
        {
            //this_step_val = end_val - begin_val
            tmp_Sk[m] = Sk[m] - tmp_Sk[m];
            tmp_Sb[m] = Sb[m] - tmp_Sb[m];
            tmp_Sn[m] = Sn[m] - tmp_Sn[m];

            //loop over tensor components
            for(int a = 0; a < 3; a++)
            {
                for(int b = 0; b <= a; b++)
                {
                    S_vals[m][(a*3)+b][s] = tmp_Sk[m].get(a,b) + tmp_Sb[m].get(a,b) + tmp_Sn[m].get(a,b);
                }
            }
        }
    }//end timestep loop
    std::cout << "|" << std::endl << std::endl; //finish progress bar

    std::ofstream outfile;

    //indices of lower-triangular matrix elements
    int indices[] = {0, 4, 8, 3, 6, 7};

    outfile.open(ofname);
    //average the stress tensor and write to file
    for(int i = 0; i < n_zvals; i++)
    {
        //Divide all components by number of timesteps to average
        for(int a = 0; a < 3; a++)
        {
            for(int b = 0; b <= a; b++)
            {
                Sb[i].set(a,b,Sb[i].get(a,b)/float(step.size()));
                Sn[i].set(a,b,Sn[i].get(a,b)/float(step.size()));
                Sk[i].set(a,b,Sk[i].get(a,b)/float(step.size()));
            }
        }
        //output to file: xx, yy, zz, yx, zx, zy
        outfile << "z= " << zvals[i] << std::endl;
        for(int m = 0; m < 6; m++)
        {
            outfile << Sk[i].get(indices[m]/3,indices[m]%3) << " ";
        }
        outfile << std::endl;
        for(int m = 0; m < 6; m++)
        {
            outfile << Sb[i].get(indices[m]/3,indices[m]%3) << " ";
        }
        outfile << std::endl;
        for(int m = 0; m < 6; m++)
        {
            outfile << Sn[i].get(indices[m]/3,indices[m]%3) << " ";
        }
        outfile << std::endl;
    }
    outfile.close();
    std::cout << "Wrote stress tensor data to: " << ofname << std::endl << std::endl;

    std::cout << "Writing timestep data to file " << tfname << std::endl;
    std::ofstream tsfile;
    tsfile.open(tfname);
    for(int i = 0; i < n_zvals; i++)
    {
        tsfile << "z= " << zvals[i] << std::endl;
        // output for xx, yy, zz, yx, zx, zy
        for(int j = 0; j < 6; j++)
        {
            for(int k = 0; k < step.size(); k++)
            {
                tsfile << S_vals[i][indices[j]][k] << " ";
            }
            tsfile << std::endl;
        }
    }
    tsfile.close();

    //clean up dynamically allocated memory
    for(unsigned int i = 0; i < step.size(); i++)
    {
        for(unsigned int j = 0; j < step[i].size(); j++)
        {
            delete[] step[i][j];
        }
    }
    for(int i = 0; i < n_zvals; i++)
    {
        for(int j = 0; j < 9; j++)
        {
            delete[] S_vals[i][j];
        }
        delete[] S_vals[i];
    }
    delete[] S_vals;
    if(NVT)
    {
        delete[] boxes[0];
    }
    else
    {
        for(unsigned int i = 0; i < boxes.size(); i++)
        {
            delete[] boxes[i];
        }
    }
    delete[] rij;
    delete[] rib;
    delete[] zvals;
    delete[] Sb;
    delete[] Sk;
    delete[] Sn;
    delete[] tmp_Sb;
    delete[] tmp_Sk;
    delete[] tmp_Sn;

    return EXIT_SUCCESS;
}


#pragma once

#include <iostream>
#include <sstream>
#include <string>
#include <cmath>
#include <vector>

//=========================================================
// Algorithm Utility Functions
//=========================================================

double fold(const double r, const double box_l)
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

double mind(const double dr, const double box_l)
{
    return dr - box_l*std::floor( (dr/box_l) + 0.5);
}

//Note: min_dr stores return value
void mind_3d(const double* ri, const double* rj, const double* box_l, double* min_dr)
{
    for(int k=0; k<3; k++)
    {
        min_dr[k] = mind(ri[k]-rj[k], box_l[k]);
    }
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

void load_data(std::string filename, std::vector<int> &type, std::vector<std::vector<int> > &bonds, std::vector<std::vector<double*> > &step, std::vector<double*> &boxes, int start=0, int stop=-1, int interval=1)
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
                type.push_back(inttemp);
                bonds.push_back(std::vector<int>() );
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
                bonds.at(inttemp).push_back(inttemp2);
            }
            else if(t[0].compare("timestep") == 0)
            {
                s++;
                if((s-start)%interval == 0 && s >= start && (s <= stop || stop == -1) )
                {
                    if( ((s-start)/interval)%update_inter == 0)
                    {
                        std::cout << "Loaded " << (s-start)/interval << " steps..." << std::endl;
                    }
                    step.push_back(std::vector<double*>() );
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
                else if((s-start)%interval == 0 && s >= start && (s <= stop || stop == -1))
                {
                    boxes.push_back(box_l);
                }
                else
                {
                    delete[] box_l;
                }
            }
            else if((s-start)%interval == 0 && s >= start && (s <= stop || stop == -1))
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
                
                step.back().push_back(coord);
            }
        }
    }
    infile.close();

    if(boxes.size() == 0)
    {
        boxes.push_back(init_box);
    }

    std::cout << "Finished Loading File.\n" << std::endl;
}

//=========================================================
// Blocking Algorithm
//=========================================================

// takes an array of doubles and its length n
// returns standard deviation of array (unbiased estimator)
double stdev(double *vals, int n)
{
    double sum = 0.0;
    double sumsq = 0.0;
    for(int i=0; i<n; i++)
    {
        sum += vals[i];
        sumsq += vals[i]*vals[i];
    }
    float N = float(n);
    return std::sqrt( (sumsq/(N-1.0)) - (sum*sum/(N*(N-1.0))) );
}

// takes an array of doubles and its length n
// blocks the array in-place, so that the array
// variable passed in now contains as its elements
// the averages of neighboring pairs of values
// that were originally in the array, and hence there
// are now half as many elements
// returns the number of elements that are now meaningful
// in the array
int block(double *vals, int n)
{
    for(int i=0; i<n/2; i++)
    {
        vals[i] = (vals[2*i] + vals[(2*i)+1]) / 2.0;
    }
    return n/2;
}

// takes an array of doubles and its length n
// returns blocked error on mean (est of std of dist of mean)
// WARNING: modifies vals IN PLACE; destroys array contents
double err_on_mean(double *vals, int n)
{
    double maxsig = 0.0;
    double sig = 0.0;
    while(n > 2)
    {
        sig = stdev(vals, n)/std::sqrt(float(n) - 1.0);
        if(sig > maxsig)
        {
            maxsig = sig;
        }
        n = block(vals, n);
    }
    return maxsig;
}

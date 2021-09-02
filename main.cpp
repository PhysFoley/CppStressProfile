#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <string>
#include <sstream>
#include <algorithm>
#include <stdlib.h>

const double pi = 3.1415926535897;

class Mat3
{
private:
    double vals[9];
public:
    Mat3();
    Mat3(double* values);
    Mat3(const Mat3& mat);
    Mat3& operator=(const Mat3& mat);
    double get(int row, int col) const;
    void set(int row, int col, double val);
};

Mat3::Mat3()
{
    for(int i=0; i<9; i++)
    {
        this->vals[i] = 0.0;
    }
}

Mat3::Mat3(const Mat3& mat)
{
    for(int i=0; i<9; i++)
    {
        this->vals[i] = mat.vals[i];
    }
}

Mat3::Mat3(double* values)
{
    for(int i=0; i<9; i++)
    {
        this->vals[i] = values[i];
    }
}

Mat3& Mat3::operator=(const Mat3& mat)
{
    for(int i=0; i<9; i++)
    {
        this->vals[i] = mat.vals[i];
    }
    return *this;
}

double Mat3::get(int row, int col) const
{
    return this->vals[row*3 + col];
}

void Mat3::set(int row, int col, double val)
{
    this->vals[row*3 + col] = val;
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


//=========================================================
// Derivatives of Interaction Potentials
//=========================================================

double (*forces[7][7])(double);

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

//head-head type interactions - 00,55,66,05,06,56
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

//head-tail interactions - 01,51,61,02,03,04,52,53,54,62,63,64
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

//attractive tail-tail interactions - 11,12,13,14,22,24,33,34,44
double f11(double r)
{
	return lj(r, 1.0, std::pow(2.0,1.0/6.0), 1.0) + cos2(r, std::pow(2.0,1.0/6.0), 1.6, 1.0);
}

//modified middle-middle interaction - not ljcos2
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
    int update_inter = 100;
    if(stop != -1)
    {
        update_inter = std::ceil(float(stop - start)/10.0);
    }
    std::ifstream infile;
    infile.open(filename);
    std::string line;
    std::vector<std::string> t;

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
                        std::cout << "Loading step " << (s-start)/interval << "..." << std::endl;
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
                
                if(s%interval == 0 && s >= start && (s <= stop || stop == -1) )
				{
					boxes->push_back(box_l);
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
            }
        }
    }
    infile.close();
    std::cout << "Finished Loading File." << std::endl;
}

//head-head type interactions - 00,55,66,05,06,56
//head-tail interactions - 01,51,61,02,03,04,52,53,54,62,63,64
//attractive tail-tail interactions - 11,12,13,14,22,24,33,34,44
int main(int argc, char** argv)
{
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
	
	

    std::stringstream ss;
    std::string ifname("centered_trajectory.vtf");
    std::string ofname("stress_profile.dat");
    std::string stdfname("profile_std.dat");
    int start = 0;
    int stop = -1;
    int interval = 1;
    int n_zvals = 51;
	double thickness = 10.0;
    double kT = 1.4;
    //double box_l[] = {17.270,17.270,20.0}; //changed this to fit my new system

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
        else if(std::string(argv[i]).compare("-s") == 0)
        {
            ss << argv[i+1];
            ss >> stdfname;
            ss.clear();
        }
        else if(std::string(argv[i]).compare("-t") == 0)
        {
            ss << argv[i+1];
            ss >> thickness;
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

    load_data(ifname, &step, &type, &bonds, &boxes, start, stop, interval);
    
    std::cout << "There were " << boxes.size() << " boxes read from the trajectory file." << std::endl;
    std::cout << "There were " << step.size() << " trajectory steps read in." << std::endl;
    /*for(int i = 0; i < step.size(); i++)
    {
        double* box_l = new double[3];
        box_l[0] = box_l[1] = 12.42;
        box_l[2] = 20.0;
        boxes.push_back(box_l);
    }*/
    
    double mean_x = 0.0;
    double mean_y = 0.0;
    for(int i=0; i < boxes.size(); i++)
    {
      mean_x += boxes[i][0];
      mean_y += boxes[i][1];
    }
    mean_x = mean_x / float(boxes.size());
    mean_y = mean_y / float(boxes.size());
    std::cout << "Mean X: " << mean_x << std::endl;
    std::cout << "Mean Y: " << mean_y << std::endl;

    double space = thickness/(n_zvals-1);
    double* zvals = new double[n_zvals];
    for(int i = 0; i < n_zvals; i++)
    {
        zvals[i]=((boxes[0][2]/2.0) - (thickness/2.0)) + i*space;
        std::cout << " " << zvals[i] << " ";
    }
    std::cout << std::endl;

    std::ofstream outfile, stdfile;
    outfile.open(ofname);
    stdfile.open(stdfname);

    Mat3* Sb = new Mat3[n_zvals];
    Mat3* Sn = new Mat3[n_zvals];
    Mat3* Sk = new Mat3[n_zvals];
    
    double r, phi_p;
    double *ri, *rj, *rb, *rij, *rib, *lo_r, *hi_r;
    double temp; // for storing temporary values
    int zbin_ind_i, zbin_ind_j, zbin_ind_b;
    
    rij = new double[3];
    rib = new double[3];
    
    double Lz = boxes[0][2]; //assumed to be constant!
    double nb_cutoff = 3.0; //make sure this is larger than longest-range nb interaction range!
    
    for(int s = 0; s < step.size(); s++)
    {
        std::cout << "Working on step " << s << std::endl;
        double A = boxes[s][0]*boxes[s][1];
        for(int i = 0; i < type.size(); i++)
        {
            //std::cout << "Particle " << i << std::endl;
            ri = step[s][i];
            zbin_ind_i = int((fold(ri[2],Lz)-zvals[0]+(space/2.0))/space);
            temp = kT/(A*space);
            for(int k = 0; k < 3; k++) //only diagonal components of kinetic tensor (ideal gas pressure)
            {
                //make sure the particle is within the domain we're treating
                if(zbin_ind_i < n_zvals && zbin_ind_i >= 0)
                {
                    Sk[zbin_ind_i].set(k,k,Sk[zbin_ind_i].get(k,k)+temp);
                }
            }
            for(int j = 0; j < i; j++) //loop over pairs (without double-counting, hence j < i)
            {
                //std::cout << "    nb interaction with particle " << j << std::endl;
                rj = step[s][j];
                mind_3d(ri, rj, boxes[s], rij); //rij stores return value
                r = std::sqrt(rij[0]*rij[0] + rij[1]*rij[1] + rij[2]*rij[2]);
                //std::cout << "    Made it past r calc..." << std::endl;
                if(r > 3.0)
                { //out of interaction range
                    continue;
                }
                else if(r == 0.0)
                {
                    std::cout << "Divide by zero in non-bonded interactions!" << std::endl;
                    continue;
                }
                phi_p = forces[type[i]][type[j]](r);
                //std::cout << "    Made it past phi_p calc..." << std::endl;
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
                                temp = rij[a]*rij[b]*phi_p/(A*space*r);
                                Sn[zbin_ind_i].set(a,b,Sn[zbin_ind_i].get(a,b)-temp);
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
                                    temp = rij[a]*rij[b]*phi_p/(A*std::abs(rij[2])*r);
                                    if(m == lo_ind) { temp = temp * lo_factor; }
                                    if(m == hi_ind) { temp = temp * hi_factor; }
                                    Sn[m].set(a,b,Sn[m].get(a,b)-temp);
                                }
                            }
                        }
                    }
                }
            } //end non-bonded interactions
            for(int bo = 0; bo < bonds[i].size(); bo++) //begin bonded interactions
            {
                //std::cout << "    bonded interaction with particle " << bonds[i][bo] << std::endl;
                rb = step[s][bonds[i][bo]];
                //std::cout << "    got rb" << std::endl;
                mind_3d(ri, rb, boxes[s], rib);
                //std::cout << "    get rib" << std::endl;
                r = std::sqrt(rib[0]*rib[0] + rib[1]*rib[1] +rib[2]*rib[2]);
                //std::cout << "    Made it past r calc..." << std::endl;
                if(r == 0.0)
                {
                    std::cout << "Divide by zero in bonded interactions!" << std::endl;
                    continue;
                }
                if(std::abs(bonds[i][bo]-i) == 1)
                {
                    phi_p = fene(r, 30.0, 1.5);
                }
                else
                {
                    phi_p = bend(r, 10.0);
                }
                //std::cout << "    Made it past phi_p calc..." << std::endl;
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
                                temp = rib[a]*rib[b]*phi_p/(A*space*r);
                                Sb[zbin_ind_i].set(a,b,Sb[zbin_ind_i].get(a,b)-temp);
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
                                    temp = 0.5*rib[a]*rib[b]*phi_p/(A*std::abs(rib[2])*r);
                                    if(m == lo_ind) { temp = temp * lo_factor; }
                                    if(m == hi_ind) { temp = temp * hi_factor; }
                                    Sb[m].set(a,b,Sb[m].get(a,b)-temp);
                                }
                            }
                        }
                    }
                }
            }//end bonded interactions
        }//end particle loop
    }//end timestep loop
    
    //Divide all components by number of timesteps to average
    for(int i = 0; i < n_zvals; i++)
    {
        for(int a = 0; a < 3; a++)
        {
            for(int b = 0; b <= a; b++)
            {
                Sb[i].set(a,b,Sb[i].get(a,b)/float(step.size()));
                Sn[i].set(a,b,Sn[i].get(a,b)/float(step.size()));
                Sk[i].set(a,b,Sk[i].get(a,b)/float(step.size()));
            }
        }
        //output to file in same style as before
        std::cout << "z= " << zvals[i] << std::endl;
        outfile << "z= " << zvals[i] << std::endl;
        for(int m = 0; m < 9; m++)
        {
            std::cout << Sk[i].get(m/3,m%3) << " ";
            outfile << Sk[i].get(m/3,m%3) << " ";
        }
        std::cout << std::endl;
        outfile << std::endl;
        for(int m = 0; m < 9; m++)
        {
            std::cout << Sb[i].get(m/3,m%3) << " ";
            outfile << Sb[i].get(m/3,m%3) << " ";
        }
        std::cout << std::endl;
        outfile << std::endl;
        for(int m = 0; m < 9; m++)
        {
            std::cout << Sn[i].get(m/3,m%3) << " ";
            outfile << Sn[i].get(m/3,m%3) << " ";
        }
        std::cout << std::endl;
        outfile << std::endl;
    }
    
    outfile.close();
    stdfile.close();

    //clean up dynamically allocated memory
    for(unsigned int i = 0; i < step.size(); i++)
    {
        for(unsigned int j = 0; j < step[i].size(); j++)
        {
            delete[] step[i][j];
        }
    }
    for(unsigned int i = 0; i < boxes.size(); i++)
    {
        delete[] boxes[i];
    }
    delete[] zvals;
    delete[] Sb;
    delete[] Sk;
    delete[] Sn;

    return 0;
}

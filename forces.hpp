#pragma once

#include <cmath>

//=========================================================
// Built-in Force Functions
//=========================================================

double lj(double r, double b, double rc, double eps)
{
    if(r <= rc)
    {
        return 24.0*eps*( (2.0*std::pow(b,12.0)/std::pow(r,13.0)) - (std::pow(b,6.0)/std::pow(r,7.0)) );
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
        return -(M_PI*eps/(2.0*wc))*std::sin(M_PI*(r-rc)/wc);
    }
}

double fene(double r, double k, double rinf)
{
    return -k*r/(1.0 - std::pow(r/rinf,2.0) );
}

double harmonic(double r, double k)
{
    return -k*(r-4.0);
}

double ljcos2(double r, double b, double rc, double wc, double eps)
{
    return lj(r, b, rc, eps) + cos2(r, rc, wc, eps);
}

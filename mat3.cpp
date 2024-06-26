#include "mat3.hpp"

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

Mat3 Mat3::operator+(const Mat3& mat) const
{
    Mat3 sum;
    for(int i=0; i<9; i++)
    {
        sum.vals[i] = this->vals[i] + mat.vals[i];
    }
    return sum;
}

Mat3 Mat3::operator-(const Mat3& mat) const
{
    Mat3 diff;
    for(int i=0; i<9; i++)
    {
        diff.vals[i] = this->vals[i] - mat.vals[i];
    }
    return diff;
}

double Mat3::get(int row, int col) const
{
    return this->vals[row*3 + col];
}

void Mat3::set(int row, int col, double val)
{
    this->vals[row*3 + col] = val;
}

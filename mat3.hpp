#pragma once

class Mat3
{
private:
    double vals[9];
public:
    Mat3();
    Mat3(double* values);
    Mat3(const Mat3& mat);
    Mat3& operator=(const Mat3& mat);
    Mat3 operator+(const Mat3& mat) const;
    Mat3 operator-(const Mat3& mat) const;
    double get(int row, int col) const;
    void set(int row, int col, double val);
};

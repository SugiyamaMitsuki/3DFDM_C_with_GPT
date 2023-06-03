// Source.h
#ifndef SOURCE_H
#define SOURCE_H
#include "array1d.h"
#include <array>
#include <cmath>

class Source
{
public:
    Source(int x, int y, int z, double MomentMagnitude, double strike, double dip, double rake);

    void setSourceLocation(int x, int y, int z);
    void setSourceParameters(double MomentMagnitude, double strike, double dip, double rake);
    void setSlipVelocity(const Array1D<double>& inputTime, const Array1D<double>& inputSlipVelocity, double dt, double totalTime);
    void makeSlipVelocity(double Vm, double td, double tb, double tr, double dt, double ts, double totalTime);
    void calculateMomentTensor(double momentTensor[3][3], double scale) const;
    void getLocation(int &x, int &y, int &z) const;
    double getSlipVelocityAtTimeStep(int index) const;
    Array1D<double> getSlipVelocity() const;

    double Radians(double degrees) const;

private:
    int x, y, z;
    double MomentMagnitude, strike, dip, rake;
    Array1D<double> slipVelocity;
};

inline double Source::Radians(double degrees) const
{
    return degrees * M_PI / 180.0;
}

#endif // SOURCE_H

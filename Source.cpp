// Source.cpp
#include "Source.h"
#include <iostream>
#include "array1d.h"

Source::Source(int x, int y, int z, double MomentMagnitude, double strike, double dip, double rake)
    : x(x), y(y), z(z), MomentMagnitude(MomentMagnitude), strike(strike), dip(dip), rake(rake) {}

void Source::setSourceLocation(int x, int y, int z)
{
    this->x = x;
    this->y = y;
    this->z = z;
}

void Source::setSourceParameters(double MomentMagnitude, double strike, double dip, double rake)
{
    this->MomentMagnitude = MomentMagnitude;
    this->strike = strike;
    this->dip = dip;
    this->rake = rake;
}

// 線形補間関数
double linearInterpolation(double x, double x1, double x2, double y1, double y2)
{
    return y1 + (x - x1) * (y2 - y1) / (x2 - x1);
}

void Source::setSlipVelocity(const Array1D<double> &inputtime, const Array1D<double> &inputSlipVelocity, double dt, double totalTime)
{
    int num_steps = totalTime / dt;
    slipVelocity.resize(0, num_steps);

    double dtInput = inputtime(1) - inputtime(0);
    double totalTimeInput = dtInput * (inputtime.size());

    // inputSlipVelocity を補間する
    if (dtInput == dt && totalTimeInput == totalTime)
    {
        slipVelocity = inputSlipVelocity; // 直接代入
    }
    else
    {
        
        Array1D<double> inputtimeInterpolated;
        inputtimeInterpolated.resize(0, num_steps);

        // std::cout << "make inputtimeInterpolated" << std::endl;

        for (int i = 0; i < num_steps; ++i)
        {
            inputtimeInterpolated(i) = i * dt;
        }

        Array1D<double> inputSlipVelocityInterpolated;
        inputSlipVelocityInterpolated.resize(0, num_steps);

        for (int i = 0; i < num_steps; ++i)
        {
            double t = inputtimeInterpolated(i);

            if (t < totalTimeInput - dtInput)
            {
                int idx = static_cast<int>(t / dtInput);
                inputSlipVelocityInterpolated(i) = linearInterpolation(t, inputtime(idx), inputtime(idx + 1), inputSlipVelocity(idx), inputSlipVelocity(idx + 1));
            }
            else
            {
                inputSlipVelocityInterpolated(i) = 0.0;
            }
            // std::cout << "i = " << i << ", t = " << t << ", inputSlipVelocityInterpolated(i) = " << inputSlipVelocityInterpolated(i) << std::endl;

        }
        slipVelocity = inputSlipVelocityInterpolated; // 直接代入
    }
}

void Source::makeSlipVelocity(double Vm, double td, double tb, double tr, double ts, double dt, double totalTime)
{
    // Calculate slip velocity
    // https://www.gms.bosai.go.jp/GMS/Documents/WebHelp/sourcetime_function4.htm

    int num_steps = totalTime / dt;
    slipVelocity.resize(0, num_steps - 1);

    double eps = tb + tb / 2 * (1.0 - tb / td / 2.0) / (1.0 - tb / td);
    if (eps > tb)
        eps = tb;
    double b = -2.0 * (2.0 * Vm / td) * (1.0 - tb / td) * pow(tb - eps, 1.5);
    double c = b / sqrt(tr - eps);
    double ar = c / (ts - tr);

    for (int i = 0; i < num_steps; ++i)
    {
        double t = i * dt;

        if (t < tb)
        {
            slipVelocity(i) = 2.0 * Vm / td * t * (1.0 - 0.5 * t / td);
        }
        else if (t < tr)
        {
            slipVelocity(i) = b / sqrt(t - eps);
        }
        else if (t < ts)
        {
            slipVelocity(i) = c - ar * (t - tr);
        }
        else
        {
            slipVelocity(i) = 0.0;
        }
    }

    // 積分値が1になるようにスケーリング
    double scale = 0.0;
    for (int i = 0; i < num_steps; ++i)
    {
        scale += slipVelocity(i) * dt;
    }
    for (int i = 0; i < num_steps; ++i)
    {
        slipVelocity(i) /= scale;
    }
}

void Source::calculateMomentTensor(double momentTensor[3][3], double scale) const
{

    // https://www.gms.bosai.go.jp/GMS/faq_item/show/3

    // X = NS (North-South) direction (positive northward)
    // Y = EW (East-West) direction (positive eastward)
    // Z = UD (Up-Down) direction (positive downward)
    // strike = φ 走向
    // dip = δ 傾斜
    // rake = λ 滑り角

    double M0 = MomentMagnitude * scale;

    // strike φ dip δ rake λ

    double Mxx = -M0 * (sin(Radians(dip)) * cos(Radians(rake)) * sin(2.0 * Radians(strike)) + sin(2.0 * Radians(dip)) * sin(Radians(rake)) * sin(Radians(strike)) * sin(Radians(strike)));
    double Myy = M0 * (sin(Radians(dip)) * cos(Radians(rake)) * sin(2.0 * Radians(strike)) - sin(2.0 * Radians(dip)) * sin(Radians(rake)) * cos(Radians(strike)) * cos(Radians(strike)));
    double Mzz = M0 * (sin(2.0 * Radians(dip)) * cos(Radians(rake)));
    double Mxy = M0 * (sin(Radians(dip)) * cos(Radians(rake)) * cos(2.0 * Radians(strike)) + sin(2.0 * Radians(dip)) * sin(Radians(rake)) * sin(2.0 * Radians(strike)) / 2.0);
    double Myz = -M0 * (cos(Radians(dip)) * cos(Radians(rake)) * sin(Radians(strike)) - cos(2.0 * Radians(dip)) * sin(Radians(rake)) * cos(Radians(strike)));
    double Mxz = -M0 * (cos(Radians(dip)) * cos(Radians(rake)) * cos(Radians(strike)) + cos(2.0 * Radians(dip)) * sin(Radians(rake)) * sin(Radians(strike)));

    momentTensor[0][0] = Mxx;
    momentTensor[0][1] = Mxy;
    momentTensor[0][2] = Mxz;
    momentTensor[1][0] = Mxy;
    momentTensor[1][1] = Myy;
    momentTensor[1][2] = Myz;
    momentTensor[2][0] = Mxz;
    momentTensor[2][1] = Myz;
    momentTensor[2][2] = Mzz;
}

double Source::getSlipVelocityAtTimeStep(int index) const
{
    return slipVelocity(index);
}
Array1D<double> Source::getSlipVelocity() const
{
    return slipVelocity;
}

void Source::getLocation(int &x, int &y, int &z) const
{
    // 震源位置を取得
    x = this->x;
    y = this->y;
    z = this->z;
}
// Grid.h
#ifndef GRID_H
#define GRID_H
#include <vector>
#include <iostream>
#include <fstream>
#include <cmath>
#include "array3d.h"
#include "array2d.h"
#include "LayerPropertyReader.h"

static const double c1 = 9.0 / 8.0;
static const double c2 = -1.0 / 24.0;

class Grid
{
public:
    Grid();

    // アクセスを提供するメソッド
    Array3D<double> &getVx() { return vx; }
    Array3D<double> &getVy() { return vy; }
    Array3D<double> &getVz() { return vz; }
    Array3D<double> &getTxx() { return txx; }
    Array3D<double> &getTyy() { return tyy; }
    Array3D<double> &getTzz() { return tzz; }
    Array3D<double> &getTyz() { return tyz; }
    Array3D<double> &getTzx() { return tzx; }
    Array3D<double> &getTxy() { return txy; }

    Array3D<double> &getRho() { return rho; }
    Array3D<double> &getMu() { return mu; }
    Array3D<double> &getLambda() { return lambda; }
    Array3D<double> &getQ0() { return Q0; }
    Array3D<double> &getPwaveSpeed() { return PwaveSpeed; }
    Array3D<double> &getSwaveSpeed() { return SwaveSpeed; }


    // 関数
    void initializeModel(int nx, int ny, int nz, double dx, double dy, double dz, double lat0, double lon0, double depth0);
    void setPhysicalProperty(bool isTest,std::string inputPhysicalProperty, std::string inputLayerProperty);
    void updateVelocity(double dt, bool useForthDerivative, double f0,
                        double spatialWeightNonReflective, double temporalWeightNonReflective);
    void applyVelocityAbsorbingBoundary(double dt, int absorbingGridNumber, double absorbingDampingFactor);


    void updateStress(double dt, double f0, bool useForthDerivative);
    void applyStressAbsorbingBoundary(double dt, int absorbingGridNumber, double absorbingDampingFactor);
    
    
    void applySource(double dt, int x, int y, int z, double momentTensor[3][3]);
    void getVelocityAtReceiver(int x, int y, int z, double &ux, double &uy, double &uz) const;

    // 諸計算
    double MeanPhysicalProperty2(double a, double b) const;
    double MeanPhysicalProperty4(double a, double b, double c, double d) const;
    double CalculateDampingFactor(double f0, double dt, double Q0val) const;

    static double UseNonreflectingCondition(double _VDtDx, double _Cx, double _Ct, double _Data_MinusX_MinusT, double _Data_PlusX_MinusT, double _Data_PlusX_PlusT);
    static double UseEdgeNonreflectingCondition(double _VDtDx, double _VDtDy, double _Cx, double _Ct,
                                                double _Data_MinusX_MinuxY_MinusT, double _Data_MinusX_PlusY_MinusT,
                                                double _Data_PlusX_MinuxY_MinusT, double _Data_MinusX_PlusY_PlusT,
                                                double _Data_PlusX__MinusY_PlusT);
    // static double WeightedData ( int _XWeight, int _YWeight, const Array3d<DataType> &_Data, int _XIndex, int _YIndex );

    // Difference operator
    double DxTxx(int i, int j, int k, bool FourthOrder) const;
    double DyTxx(int i, int j, int k, bool FourthOrder) const;
    double DzTxx(int i, int j, int k, bool FourthOrder) const;
    double DxTyy(int i, int j, int k, bool FourthOrder) const;
    double DyTyy(int i, int j, int k, bool FourthOrder) const;
    double DzTyy(int i, int j, int k, bool FourthOrder) const;
    double DxTzz(int i, int j, int k, bool FourthOrder) const;
    double DyTzz(int i, int j, int k, bool FourthOrder) const;
    double DzTzz(int i, int j, int k, bool FourthOrder) const;

    double DyTyz(int i, int j, int k, bool FourthOrder) const;
    double DzTyz(int i, int j, int k, bool FourthOrder) const;
    double DzTzx(int i, int j, int k, bool FourthOrder) const;
    double DxTzx(int i, int j, int k, bool FourthOrder) const;
    double DxTxy(int i, int j, int k, bool FourthOrder) const;
    double DyTxy(int i, int j, int k, bool FourthOrder) const;

    double DxVx(int i, int j, int k, bool FourthOrder) const;
    double DxVy(int i, int j, int k, bool FourthOrder) const;
    double DxVz(int i, int j, int k, bool FourthOrder) const;
    double DyVx(int i, int j, int k, bool FourthOrder) const;
    double DyVy(int i, int j, int k, bool FourthOrder) const;
    double DyVz(int i, int j, int k, bool FourthOrder) const;
    double DzVx(int i, int j, int k, bool FourthOrder) const;
    double DzVy(int i, int j, int k, bool FourthOrder) const;
    double DzVz(int i, int j, int k, bool FourthOrder) const;

    // Boundary condition
    // double TxxAtBoundary(int i, int j, int k) const;

private:
    int nx, ny, nz;
    double dx, dy, dz;
    double lat0, lon0, depth0;

    // 3D array
    Array3D<double> lambda;
    Array3D<double> mu;
    Array3D<double> rho;
    Array3D<double> Q0;
    Array3D<double> PwaveSpeed;
    Array3D<double> SwaveSpeed;

    Array3D<double> vx;
    Array3D<double> vy;
    Array3D<double> vz;

    Array3D<double> txx;
    Array3D<double> tyy;
    Array3D<double> tzz;

    Array3D<double> tyz;
    Array3D<double> tzx;
    Array3D<double> txy;

    // 境界層の配列
    Array2D<double> bufferVx_MinusX;
    Array2D<double> bufferVy_MinusX;
    Array2D<double> bufferVz_MinusX;
    Array2D<double> bufferVx_PlusX;
    Array2D<double> bufferVy_PlusX;
    Array2D<double> bufferVz_PlusX;

    Array2D<double> bufferVx_MinusY;
    Array2D<double> bufferVy_MinusY;
    Array2D<double> bufferVz_MinusY;
    Array2D<double> bufferVx_PlusY;
    Array2D<double> bufferVy_PlusY;
    Array2D<double> bufferVz_PlusY;

    Array2D<double> bufferVx_PlusZ;
    Array2D<double> bufferVy_PlusZ;
    Array2D<double> bufferVz_PlusZ;
};

// mean physical property
inline double Grid::MeanPhysicalProperty2(double a, double b) const
{
    return (a + b) / 2.0;
}
inline double Grid::MeanPhysicalProperty4(double a, double b, double c, double d) const
{
    return (a + b + c + d) / 4.0;
}
// calculate damping factor
inline double Grid::CalculateDampingFactor(double f0, double dt, double Q0val) const
{   
    return exp(-M_PI * f0 * dt / Q0val);
}

// 無反射境界条件
inline double Grid::UseNonreflectingCondition(double vDtDx, double Cx, double Ct, double data_MinusX_MinusT, double data_PlusX_MinusT, double data_PlusX_PlusT)
{
    double a11 = Cx - vDtDx * Ct;
    double a12 = Cx + vDtDx * (1.0 - Ct);
    double a21 = (1.0 - Cx) + vDtDx * Ct;
    double a22 = (1.0 - Cx) - vDtDx * (1.0 - Ct);

    return (a11 * data_MinusX_MinusT + a21 * data_PlusX_MinusT - a22 * data_PlusX_PlusT) / a12;
}
// 無反射境界条件
inline double Grid::UseEdgeNonreflectingCondition(double vDtDx, double vDtDy, double Cx, double Ct, double data_MinusX_MinusY_MinusT, double data_MinusX_PlusY_MinusT, double data_PlusX__MinusY_MinusT, double data_MinusX_PlusY_PlusT, double data_PlusX_MinusY_PlusT)
{
    double a111 = 2.0 * Cx - (vDtDx + vDtDy) * Ct;
    double a112 = 2.0 * Cx - (vDtDx + vDtDy) * (1.0 - Ct);
    double a121 = (1.0 - Cx) + vDtDx * Ct;
    double a211 = (1.0 - Cx) + vDtDy * Ct;
    double a122 = (1.0 - Cx) - vDtDx * (1.0 - Ct);
    double a212 = (1.0 - Cx) - vDtDy * (1.0 - Ct);

    return (a111 * data_MinusX_MinusY_MinusT + a121 * data_MinusX_PlusY_MinusT + a211 * data_PlusX__MinusY_MinusT - a122 * data_MinusX_PlusY_PlusT - a212 * data_PlusX_MinusY_PlusT) / a112;
}

// Difference operator
// Txx, Tyy, Tzz
inline double Grid::DxTxx(int i, int j, int k, bool FourthOrder) const
{
    if (FourthOrder && (0 < i && i < nx - 2))
    {
        return (c1 * (txx(i + 1, j, k) - txx(i, j, k)) + c2 * (txx(i + 2, j, k) - txx(i - 1, j, k))) / dx;
    }
    else
    {
        return (txx(i + 1, j, k) - txx(i, j, k)) / dx;
    }
}
inline double Grid::DyTxx(int i, int j, int k, bool FourthOrder) const
{
    if (FourthOrder && (0 < j && j < ny - 2))
    {
        return (c1 * (txx(i, j + 1, k) - txx(i, j, k)) + c2 * (txx(i, j + 2, k) - txx(i, j - 1, k))) / dy;
    }
    else
    {
        return (txx(i, j + 1, k) - txx(i, j, k)) / dy;
    }
}
inline double Grid::DzTxx(int i, int j, int k, bool FourthOrder) const
{
    if (FourthOrder && (0 < k && k < nz - 2))
    {
        return (c1 * (txx(i, j, k + 1) - txx(i, j, k)) + c2 * (txx(i, j, k + 2) - txx(i, j, k - 1))) / dz;
    }
    else
    {
        return (txx(i, j, k + 1) - txx(i, j, k)) / dz;
    }
}
inline double Grid::DxTyy(int i, int j, int k, bool FourthOrder) const
{
    if (FourthOrder && (0 < i && i < nx - 2))
    {
        return (c1 * (tyy(i + 1, j, k) - tyy(i, j, k)) + c2 * (tyy(i + 2, j, k) - tyy(i - 1, j, k))) / dx;
    }
    else
    {
        return (tyy(i + 1, j, k) - tyy(i, j, k)) / dx;
    }
}
inline double Grid::DyTyy(int i, int j, int k, bool FourthOrder) const
{
    if (FourthOrder && (0 < j && j < ny - 2))
    {
        return (c1 * (tyy(i, j + 1, k) - tyy(i, j, k)) + c2 * (tyy(i, j + 2, k) - tyy(i, j - 1, k))) / dy;
    }
    else
    {
        return (tyy(i, j + 1, k) - tyy(i, j, k)) / dy;
    }
}
inline double Grid::DzTyy(int i, int j, int k, bool FourthOrder) const
{
    if (FourthOrder && (0 < k && k < nz - 2))
    {
        return (c1 * (tyy(i, j, k + 1) - tyy(i, j, k)) + c2 * (tyy(i, j, k + 2) - tyy(i, j, k - 1))) / dz;
    }
    else
    {
        return (tyy(i, j, k + 1) - tyy(i, j, k)) / dz;
    }
}
inline double Grid::DxTzz(int i, int j, int k, bool FourthOrder) const
{
    if (FourthOrder && (0 < i && i < nx - 2))
    {
        return (c1 * (tzz(i + 1, j, k) - tzz(i, j, k)) + c2 * (tzz(i + 2, j, k) - tzz(i - 1, j, k))) / dx;
    }
    else
    {
        return (tzz(i + 1, j, k) - tzz(i, j, k)) / dx;
    }
}
inline double Grid::DyTzz(int i, int j, int k, bool FourthOrder) const
{
    if (FourthOrder && (0 < j && j < ny - 2))
    {
        return (c1 * (tzz(i, j + 1, k) - tzz(i, j, k)) + c2 * (tzz(i, j + 2, k) - tzz(i, j - 1, k))) / dy;
    }
    else
    {
        return (tzz(i, j + 1, k) - tzz(i, j, k)) / dy;
    }
}
inline double Grid::DzTzz(int i, int j, int k, bool FourthOrder) const
{
    if (FourthOrder && (0 < k && k < nz - 2))
    {
        return (c1 * (tzz(i, j, k + 1) - tzz(i, j, k)) + c2 * (tzz(i, j, k + 2) - tzz(i, j, k - 1))) / dz;
    }
    else
    {
        return (tzz(i, j, k + 1) - tzz(i, j, k)) / dz;
    }
}


// Txy, Tyx, Txz, Tzx, Tyz, Tzy
inline double Grid::DyTyz(int i, int j, int k, bool FourthOrder) const
{
    if (FourthOrder && (1 < j && j < ny - 1))
    {
        return (c1 * (tyz(i, j, k) - tyz(i, j - 1, k)) + c2 * (tyz(i, j + 1, k) - tyz(i, j - 2, k))) / dy;
    }
    else
    {
        return (tyz(i, j, k) - tyz(i, j - 1, k)) / dy;
    }
}
inline double Grid::DzTyz(int i, int j, int k, bool FourthOrder) const
{
    if (FourthOrder && (1 < k && k < nz - 1))
    {
        return (c1 * (tyz(i, j, k) - tyz(i, j, k - 1)) + c2 * (tyz(i, j, k + 1) - tyz(i, j, k - 2))) / dz;
    }
    else
    {
        return (tyz(i, j, k) - tyz(i, j, k - 1)) / dz;
    }
}
inline double Grid::DzTzx(int i, int j, int k, bool FourthOrder) const
{
    if (FourthOrder && (1 < k && k < nz - 1))
    {
        return (c1 * (tzx(i, j, k) - tzx(i, j, k - 1)) + c2 * (tzx(i, j, k + 1) - tzx(i, j, k - 2))) / dz;
    }
    else
    {
        return (tzx(i, j, k) - tzx(i, j, k - 1)) / dz;
    }
}
inline double Grid::DxTzx(int i, int j, int k, bool FourthOrder) const
{
    if (FourthOrder && (1 < i && i < nx - 1))
    {
        return (c1 * (tzx(i, j, k) - tzx(i - 1, j, k)) + c2 * (tzx(i + 1, j, k) - tzx(i - 2, j, k))) / dx;
    }
    else
    {
        return (tzx(i, j, k) - tzx(i - 1, j, k)) / dx;
    }
}
inline double Grid::DxTxy(int i, int j, int k, bool FourthOrder) const
{
    if (FourthOrder && (1 < i && i < nx - 1))
    {
        return (c1 * (txy(i, j, k) - txy(i - 1, j, k)) + c2 * (txy(i + 1, j, k) - txy(i - 2, j, k))) / dx;
    }
    else
    {
        return (txy(i, j, k) - txy(i - 1, j, k)) / dx;
    }
}
inline double Grid::DyTxy(int i, int j, int k, bool FourthOrder) const
{
    if (FourthOrder && (1 < j && j < ny - 1))
    {
        return (c1 * (txy(i, j, k) - txy(i, j - 1, k)) + c2 * (txy(i, j + 1, k) - txy(i, j - 2, k))) / dy;
    }
    else
    {
        return (txy(i, j, k) - txy(i, j - 1, k)) / dy;
    }
}


// Difference operator
// dxvx dyvy dzvz
inline double Grid::DxVx(int i, int j, int k, bool FourthOrder) const
{
    if (FourthOrder && (1 < i && i < nx - 1))
    {
        return (c1 * (vx(i, j, k) - vx(i - 1, j, k)) + c2 * (vx(i + 1, j, k) - vx(i - 2, j, k))) / dx;
    }
    else
    {
        return (vx(i, j, k) - vx(i - 1, j, k)) / dx;
    }
}
inline double Grid::DyVy(int i, int j, int k, bool FourthOrder) const
{
    if (FourthOrder && (1 < j && j < ny - 1))
    {
        return (c1 * (vy(i, j, k) - vy(i, j - 1, k)) + c2 * (vy(i, j + 1, k) - vy(i, j - 2, k))) / dy;
    }
    else
    {
        return (vy(i, j, k) - vy(i, j - 1, k)) / dy;
    }
}
inline double Grid::DzVz(int i, int j, int k, bool FourthOrder) const
{
    if (FourthOrder && (1 < k && k < nz - 1))
    {
        return (c1 * (vz(i, j, k) - vz(i, j, k - 1)) + c2 * (vz(i, j, k + 1) - vz(i, j, k - 2))) / dz;
    }
    else
    {
        return (vz(i, j, k) - vz(i, j, k - 1)) / dz;
    }
}

// dxvy dxvz dyvx dyvz dzvx dzvy
inline double Grid::DxVy(int i, int j, int k, bool FourthOrder) const
{
    if (FourthOrder && (0 < i && i < nx - 2))
    {
        return (c1 * (vy(i + 1, j, k) - vy(i, j, k)) + c2 * (vy(i + 2, j, k) - vy(i - 1, j, k))) / dx;
    }
    else
    {
        return (vy(i + 1, j, k) - vy(i, j, k)) / dx;
    }
}
inline double Grid::DxVz(int i, int j, int k, bool FourthOrder) const
{
    if (FourthOrder && (0 < i && i < nx - 2))
    {
        return (c1 * (vz(i + 1, j, k) - vz(i, j, k)) + c2 * (vz(i + 2, j, k) - vz(i - 1, j, k))) / dx;
    }
    else
    {
        return (vz(i + 1, j, k) - vz(i, j, k)) / dx;
    }
}
inline double Grid::DyVx(int i, int j, int k, bool FourthOrder) const
{
    if (FourthOrder && (0 < j && j < ny - 2))
    {
        return (c1 * (vx(i, j + 1, k) - vx(i, j, k)) + c2 * (vx(i, j + 2, k) - vx(i, j - 1, k))) / dy;
    }
    else
    {
        return (vx(i, j + 1, k) - vx(i, j, k)) / dy;
    }
}

inline double Grid::DyVz(int i, int j, int k, bool FourthOrder) const
{
    if (FourthOrder && (0 < j && j < ny - 2))
    {
        return (c1 * (vz(i, j + 1, k) - vz(i, j, k)) + c2 * (vz(i, j + 2, k) - vz(i, j - 1, k))) / dy;
    }
    else
    {
        return (vz(i, j + 1, k) - vz(i, j, k)) / dy;
    }
}
inline double Grid::DzVx(int i, int j, int k, bool FourthOrder) const
{
    if (FourthOrder && (0 < k && k < nz - 2))
    {
        return (c1 * (vx(i, j, k + 1) - vx(i, j, k)) + c2 * (vx(i, j, k + 2) - vx(i, j, k - 1))) / dz;
    }
    else
    {
        return (vx(i, j, k + 1) - vx(i, j, k)) / dz;
    }
}
inline double Grid::DzVy(int i, int j, int k, bool FourthOrder) const
{
    if (FourthOrder && (0 < k && k < nz - 2))
    {
        return (c1 * (vy(i, j, k + 1) - vy(i, j, k)) + c2 * (vy(i, j, k + 2) - vy(i, j, k - 1))) / dz;
    }
    else
    {
        return (vy(i, j, k + 1) - vy(i, j, k)) / dz;
    }
}

#endif // GRID_H

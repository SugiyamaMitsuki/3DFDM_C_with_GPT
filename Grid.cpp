// Grid.cpp
#include "Grid.h"
#include <cmath>
#include <iostream>
#include <fstream>
#include "array1d.h"
#include "MaterialPropertyReader.h"
#include "CoordinateManager.h"
#include "LayerPropertyReader.h"

Grid::Grid() {}

// memo
// 物性値のデータ数は、格子点の数と同じ
// 速度のデータ数は、成分方向が格子点の数より1多い
// 軸応力のデータ数は、格子点の数と同じ
// せん断応力のデータ数は、成分方向が格子点の数より1多い

void Grid::initializeModel(int nx, int ny, int nz, double dx, double dy, double dz, double lat0, double lon0, double depth0)
{
    this->nx = nx;
    this->ny = ny;
    this->nz = nz;
    this->dx = dx;
    this->dy = dy;
    this->dz = dz;
    this->lat0 = lat0;
    this->lon0 = lon0;
    this->depth0 = depth0;

    std::cout << "nx: " << nx << std::endl;
    std::cout << "ny: " << ny << std::endl;
    std::cout << "nz: " << nz << std::endl;
    std::cout << "dx: " << dx << std::endl;
    std::cout << "dy: " << dy << std::endl;
    std::cout << "dz: " << dz << std::endl;
    std::cout << "lat0: " << lat0 << std::endl;
    std::cout << "lon0: " << lon0 << std::endl;
    std::cout << "depth0: " << depth0 << std::endl;

    // lambda(nx, ny, nz);
    // mu(nx, ny, nz);
    // rho(nx, ny, nz);
    // Q0(nx, ny, nz);
    // PwaveSpeed(nx, ny, nz);
    // SwaveSpeed(nx, ny, nz);
    // vx(nx + 1, ny, nz);
    // vy(nx, ny + 1, nz);
    // vz(nx, ny, nz + 1);
    // txx(nx, ny, nz);
    // tyy(nx, ny, nz);
    // tzz(nx, ny, nz);
    // tyz(nx + 1, ny, nz);
    // tzx(nx, ny + 1, nz);
    // txy(nx, ny, nz + 1);
    // bufferVx_MinusX(ny, nz);
    // bufferVy_MinusX(ny - 1, nz);
    // bufferVz_MinusX(ny, nz - 1);
    // bufferVx_PlusX(ny, nz);
    // bufferVy_PlusX(ny - 1, nz);
    // bufferVz_PlusX(ny, nz - 1);
    // bufferVx_MinusY(nz, nx - 1);
    // bufferVy_MinusY(nz, nx);
    // bufferVz_MinusY(nz - 1, nx);
    // bufferVx_PlusY(nz, nx - 1);
    // bufferVy_PlusY(nz, nx);
    // bufferVz_PlusY(nz - 1, nx);
    // bufferVx_PlusZ(nx + 1, ny);
    // bufferVy_PlusZ(nx, ny + 1);
    // bufferVz_PlusZ(nx, ny);

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // インデックスの調整
    // 食い違い格子に対応する
    //      resizeはアクセスしたいインデックスの範囲を指定する（0, nx - 1）--> nx個の格子点
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // physical property
    lambda.resize(0, nx - 1, 0, ny - 1, 0, nz - 1);
    mu.resize(0, nx - 1, 0, ny - 1, 0, nz - 1);
    rho.resize(0, nx - 1, 0, ny - 1, 0, nz - 1);
    Q0.resize(0, nx - 1, 0, ny - 1, 0, nz - 1);
    PwaveSpeed.resize(0, nx - 1, 0, ny - 1, 0, nz - 1);
    SwaveSpeed.resize(0, nx - 1, 0, ny - 1, 0, nz - 1);

    // velocity 成分方向が格子点の数より1多い
    vx.resize(-1, nx - 1, 0, ny - 1, 0, nz - 1);
    vy.resize(0, nx - 1, -1, ny - 1, 0, nz - 1);
    vz.resize(0, nx - 1, 0, ny - 1, -1, nz - 1);

    // stress　軸応力のデータ数は、格子点の数と同じ
    txx.resize(0, nx - 1, 0, ny - 1, 0, nz - 1);
    tyy.resize(0, nx - 1, 0, ny - 1, 0, nz - 1);
    tzz.resize(0, nx - 1, 0, ny - 1, 0, nz - 1);

    // shear stress　せん断応力のデータ数は、成分方向が格子点の数より1多い
    tyz.resize(-1, nx - 1, 0, ny - 1, 0, nz - 1);
    tzx.resize(0, nx - 1, -1, ny - 1, 0, nz - 1);
    txy.resize(0, nx - 1, 0, ny - 1, -1, nz - 1);

    // buffer
    bufferVx_MinusX.resize(0, ny - 1, 0, nz - 1);
    bufferVy_MinusX.resize(0, ny - 2, 0, nz - 1);
    bufferVz_MinusX.resize(0, ny - 1, 0, nz - 2);
    bufferVx_PlusX.resize(0, ny - 1, 0, nz - 1);
    bufferVy_PlusX.resize(0, ny - 2, 0, nz - 1);
    bufferVz_PlusX.resize(0, ny - 1, 0, nz - 2);

    bufferVx_MinusY.resize(0, nz - 1, 0, nx - 2);
    bufferVy_MinusY.resize(0, nz - 1, 0, nx - 1);
    bufferVz_MinusY.resize(0, nz - 2, 0, nx - 1);
    bufferVx_PlusY.resize(0, nz - 1, 0, nx - 2);
    bufferVy_PlusY.resize(0, nz - 1, 0, nx - 1);
    bufferVz_PlusY.resize(0, nz - 2, 0, nx - 1);

    bufferVx_PlusZ.resize(-1, nx - 1, 0, ny - 1);
    bufferVy_PlusZ.resize(0, nx - 1, -1, ny - 1);
    bufferVz_PlusZ.resize(0, nx - 1, 0, ny - 1);
}

void Grid::setPhysicalProperty(bool isTest, std::string inputPhysicalProperty, std::string inputLayerProperty)
{
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // 地盤モデルを作成
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // 地盤構造の設定: 物性値 (lambda, mu, rho) を設定

    // Test運用と本番運用で切り替える

    if (isTest) // Test 成層構造
    {

        // Test
        // ここでは、3層の地盤モデルを作成します。各層の厚さ、Vp、Vs、および密度を設定します。
        std::vector<double> layerThickness = {2500.0, 2000.0, 250000.0}; // m
        std::vector<double> layerVp = {3400.0, 3600.0, 3800.0};          // m/s
        std::vector<double> layerVs = {1600.0, 1800.0, 2000.0};          // m/s
        std::vector<double> layerDensity = {2300.0, 2350.0, 2400.0};     // kg/m^3
        std::vector<double> layerQ0 = {150.0, 150.0, 200.0};             // Q0

        // 各グリッド点の物性値を設定します。
        for (int i = 0; i < nx; ++i)
        {
            for (int j = 0; j < ny; ++j)
            {
                double depth = 0.0;
                for (int k = 0; k < nz; ++k)
                {
                    // debug
                    // std::cout << "i = " << i << ", j = " << j << ", k = " << k << std::endl;

                    // 現在の深さに対応する層を見つけます。
                    int layerIndex = 0;
                    double cumulativeThickness = 0.0;
                    while (cumulativeThickness < depth && layerIndex + 1 < layerThickness.size())
                    {
                        layerIndex++;
                        cumulativeThickness += layerThickness[layerIndex];
                    }

                    // Set the material properties.
                    double vp = layerVp[layerIndex];
                    double vs = layerVs[layerIndex];
                    double rhoVal = layerDensity[layerIndex];
                    double muVal = rhoVal * vs * vs;
                    double lambdaVal = rhoVal * vp * vp - 2.0 * muVal;
                    double Q0val = layerQ0[layerIndex];

                    // lambda, mu, rhoを設定します。
                    lambda(i, j, k) = lambdaVal;
                    mu(i, j, k) = muVal;
                    rho(i, j, k) = rhoVal;
                    // Q0を設定します.
                    Q0(i, j, k) = Q0val;
                    // 速度を設定します。
                    PwaveSpeed(i, j, k) = vp;
                    SwaveSpeed(i, j, k) = vs;

                    // 次の深さを計算します。
                    depth += dz;
                }
            }
        }
    }
    else // 本番
    {
        // J-SHS 地盤モデルを使用

        // debug
        // std::cout << "inputPhysicalProperty = " << inputPhysicalProperty << std::endl;
        // std::cout << "inputLayerProperty = " << inputLayerProperty << std::endl;
        // std::cout << "nx = " << nx << ", ny = " << ny << ", nz = " << nz << std::endl;
        // std::cout << "dx = " << dx << ", dy = " << dy << ", dz = " << dz << std::endl;

        // 物性値の読み込み
        MaterialPropertyReader materialproperty(inputPhysicalProperty);
        // MaterialProperty property = reader[stn];
        // STN: 物性値番号
        // SVP: P 波速度(m/s)
        // SVS: S 波速度(m/s)
        // SRO: 密度(kg/m3)
        // SQP: Qp 値※ (1Hz における値)
        // SQS: Qs 値※ (1Hz における値)

        // Layerデータの読み込み
        // X方向長さ (m)
        double xLength = nx * dx;
        // Y方向長さ (m)
        double yLength = ny * dy;
        // Z方向長さ (m)
        double zLength = nz * dz;

        // 緯度経度の範囲を取得します。
        double deltaLat = CoordinateManager::GetLatDifference(xLength, lat0);
        double deltaLon = CoordinateManager::GetLonDifference(yLength, lat0);
        // maxLat, minLat, maxLon, minLonを計算
        double maxLat = lat0 + deltaLat / 2.0;
        double minLat = lat0 - deltaLat / 2.0;
        double maxLon = lon0 + deltaLon / 2.0;
        double minLon = lon0 - deltaLon / 2.0;

        // Layerデータを読み込みます
        double margin = 0.1;
        LayerPropertyReader layerPropertyReader(inputLayerProperty, minLat - margin, maxLat + margin, minLon - margin, maxLon + margin);

        // 各グリッド点の物性値を設定します。
        for (int i = 0; i < nx; ++i)
        {
            for (int j = 0; j < ny; ++j)
            {
                double targetLat = lat0 - (maxLat - minLat) / 2.0 + static_cast<double>(i) / nx * deltaLat;
                double targetLon = lon0 - (maxLon - minLon) / 2.0 + static_cast<double>(j) / ny * deltaLon;

                // i,jに対応するメッシュコードを計算します。
                std::string targetMeshCode = CoordinateManager::LatLonToMeshCode(targetLat, targetLon);

                // メッシュコードに対応するレイヤー情報（LayerProperty）を取得します。
                LayerProperty targetLayerProperty = layerPropertyReader[targetMeshCode.substr(0, 8)];

                double depth = depth0; // 基準深さを起点として、深さを計算します。
                for (int k = 0; k < nz; ++k)
                {
                    // レイヤー情報（LayerProperty）から、深さに対応する物性値を取得します。
                    double elevation = depth * (-1.0); // 標高　上向きが正
                    int targetSTN = layerPropertyReader.getLayerByElevation(targetLayerProperty, elevation);

                    // Set the material properties.
                    MaterialProperty targetLayer = materialproperty[targetSTN];

                    double vp = targetLayer.SVP;
                    double vs = targetLayer.SVS;
                    double rhoVal = targetLayer.SRO;
                    double muVal = rhoVal * vs * vs;
                    double lambdaVal = rhoVal * vp * vp - 2.0 * muVal;
                    double Q0val = targetLayer.SQS;

                    // lambda, mu, rhoを設定します。
                    lambda(i, j, k) = lambdaVal;
                    mu(i, j, k) = muVal;
                    rho(i, j, k) = rhoVal;
                    // Q0を設定します.
                    Q0(i, j, k) = Q0val;
                    // 速度を設定します。
                    PwaveSpeed(i, j, k) = vp;
                    SwaveSpeed(i, j, k) = vs;

                    // 次の深さを計算します。
                    depth += dz;


                    // debug
                    // std::cout << "i = " << i << ", j = " << j << ", k = " << k << std::endl;
                    // std::cout << "targetLat = " << targetLat << ", targetLon = " << targetLon << std::endl;
                    // std::cout << "targetMeshCode = " << targetMeshCode << std::endl;
                    // std::cout << "targetSTN = " << targetSTN << std::endl;
                    // std::cout << "elevation = " << elevation << std::endl;
                    // std::cout << "vp = " << vp << ", vs = " << vs << ", rhoVal = " << rhoVal << std::endl;
                    // std::cout << "muVal = " << muVal << ", lambdaVal = " << lambdaVal << std::endl;
                    // std::cout << "Q0val = " << Q0val << std::endl;
                    
                }
            }
        }
    }
}

void Grid::updateVelocity(double dt, bool useForthDerivative,
                          double f0,
                          double spatialWeightNonReflective, double temporalWeightNonReflective)
{

    // 速度更新式に基づいて、各グリッド点の速度成分を計算
    // ... (Graves, 1996およびPitarka, 1999の速度更新式を参照)

    // std::cout << "updateVelocity" << std::endl;
    // std::cout << "" << std::endl;

#pragma omp parallel
    {
        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        //
        // 境界データ 境界の一つ内側の値をbufferに格納する / 平井先生のコードを参考にしています。
        //
        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

        // YZ平面
#pragma omp for collapse(2)
        for (int j = 0; j < ny; ++j)
        {
            for (int k = 0; k < nz - 1; ++k)
            {

                bufferVx_MinusX(j, k) = vx(0, j, k);
                bufferVx_PlusX(j, k) = vx(nx - 2, j, k);
            }
        }
#pragma omp for collapse(2)
        for (int j = 0; j < ny - 1; ++j)
        {
            for (int k = 0; k < nz - 1; ++k)
            {
                bufferVy_MinusX(j, k) = vy(1, j, k);
                bufferVy_PlusX(j, k) = vy(nx - 2, j, k);
            }
        }
#pragma omp for collapse(2)
        for (int j = 0; j < ny; ++j)
        {
            for (int k = 0; k < nz - 1; ++k)
            {
                bufferVz_MinusX(j, k) = vz(1, j, k);
                bufferVz_PlusX(j, k) = vz(nx - 2, j, k);
            }
        }

        // XZ平面
#pragma omp for collapse(2)
        for (int k = 0; k < nz - 1; ++k)
        {
            for (int i = 0; i < nx - 1; ++i)
            {
                bufferVx_MinusY(k, i) = vx(i, 1, k);
                bufferVx_PlusY(k, i) = vx(i, ny - 2, k);
            }
        }
#pragma omp for collapse(2)
        for (int k = 0; k < nz - 1; ++k)
        {
            for (int i = 0; i < nx; ++i)
            {
                bufferVy_MinusY(k, i) = vy(i, 0, k);
                bufferVy_PlusY(k, i) = vy(i, ny - 2, k);
            }
        }
#pragma omp for collapse(2)
        for (int k = 0; k < nz - 1; ++k)
        {
            for (int i = 0; i < nx; ++i)
            {
                bufferVz_MinusY(k, i) = vz(i, 1, k);
                bufferVz_PlusY(k, i) = vz(i, ny - 2, k);
            }
        }

        // XY平面
#pragma omp for collapse(2)
        for (int i = -1; i < nx; ++i)
        {
            for (int j = 0; j < ny; ++j)
            {
                bufferVx_PlusZ(i, j) = vx(i, j, nz - 2);
            }
        }
#pragma omp for collapse(2)
        for (int i = 0; i < nx; ++i)
        {
            for (int j = -1; j < ny; ++j)
            {
                bufferVy_PlusZ(i, j) = vy(i, j, nz - 2);
            }
        }
#pragma omp for collapse(2)
        for (int i = 0; i < nx; ++i)
        {
            for (int j = 0; j < ny; ++j)
            {
                bufferVz_PlusZ(i, j) = vz(i, j, nz - 2);
            }
        }

        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        //
        // 境界面以外の計算　境界から一層内側は2次精度、それ以外は4次精度
        //
        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        //
        // x方向の速度更新式 vx(-1, nx,   0, ny,   0, nz)
        //
        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#pragma omp for collapse(3)
        for (int i = 0; i < nx - 1; ++i)
        {
            for (int j = 1; j < ny - 1; ++j)
            {
                for (int k = 0; k < nz - 1; ++k)
                {
                    double buoyancy = 1.0 / Grid::MeanPhysicalProperty2(rho(i, j, k), rho(i + 1, j, k));
                    double Q0val = Grid::MeanPhysicalProperty2(Q0(i, j, k), Q0(i + 1, j, k));
                    double damping = Grid::CalculateDampingFactor(f0, dt, Q0val);

                    double dxTxx = DxTxx(i, j, k, useForthDerivative);
                    double dyTxy = DyTxy(i, j, k, useForthDerivative);
                    double dzTxz = (k == 0) ? 2.0 * tzx(i, j, 0) / dz : DzTzx(i, j, k, useForthDerivative);

                    vx(i, j, k) = damping * (vx(i, j, k) + buoyancy * (dxTxx + dyTxy + dzTxz) * dt);
                }
            }
        }
        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        //
        // y方向の速度更新式 vy
        //
        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#pragma omp for collapse(3)
        for (int i = 1; i < nx - 1; ++i)
        {
            for (int j = 0; j < ny - 1; ++j)
            {
                for (int k = 0; k < nz - 1; ++k)
                {
                    double buoyancy = 1.0 / Grid::MeanPhysicalProperty2(rho(i, j, k), rho(i, j + 1, k));
                    double Q0val = Grid::MeanPhysicalProperty2(Q0(i, j, k), Q0(i, j + 1, k));
                    double damping = Grid::CalculateDampingFactor(f0, dt, Q0val);

                    double dxTxy = DxTxy(i, j, k, useForthDerivative);
                    double dyTyy = DyTyy(i, j, k, useForthDerivative);
                    double dzTyz = (k == 0) ? 2.0 * tyz(i, j, 0) / dz : DzTyz(i, j, k, useForthDerivative);

                    vy(i, j, k) = damping * (vy(i, j, k) + buoyancy * (dxTxy + dyTyy + dzTyz) * dt);
                }
            }
        }
        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        //
        // z方向の速度更新式 vz
        //
        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#pragma omp for collapse(3)
        for (int i = 1; i < nx - 1; ++i)
        {
            for (int j = 1; j < ny - 1; ++j)
            {
                for (int k = 0; k < nz - 1; ++k)
                {
                    double buoyancy = 1.0 / Grid::MeanPhysicalProperty2(rho(i, j, k), rho(i, j, k + 1));
                    double Q0val = Grid::MeanPhysicalProperty2(Q0(i, j, k), Q0(i, j, k + 1));
                    double damping = Grid::CalculateDampingFactor(f0, dt, Q0val);

                    double dxTxz = DxTzx(i, j, k, useForthDerivative);
                    double dyTyz = DyTyz(i, j, k, useForthDerivative);
                    double dzTzz = DzTzz(i, j, k, useForthDerivative);

                    vz(i, j, k) = damping * (vz(i, j, k) + buoyancy * (dxTxz + dyTyz + dzTzz) * dt);
                }
            }
        }

        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        //
        // 無反射境界条件 （地表面以外の境界）
        //
        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

        // 　1.Vy,VzをX方向に拡大
        // 　2.Vx,VzをY方向に拡大
        // 　3.Vx,VyをZ正方向に拡大
        // 　4.VzをXY平面の四隅へ拡大
        // 　5.VxをX方向へ拡大
        // 　6.VyをY方向へ拡大
        // 　7.Vx,Vy,VzをZ正負方向に拡大

        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        // 1.Vy,VzをX方向に拡大
        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#pragma omp for collapse(2)
        for (int j = 0; j < ny - 1; ++j)
        {
            for (int k = 0; k < nz - 1; ++k)
            {
                double vs = Grid::MeanPhysicalProperty4(SwaveSpeed(0, j, k), SwaveSpeed(0, j + 1, k), SwaveSpeed(1, j, k), SwaveSpeed(1, j + 1, k));
                double vDtDx = vs * dt / dx;
                vy(0, j, k) = Grid::UseNonreflectingCondition(vDtDx, spatialWeightNonReflective, temporalWeightNonReflective, vy(0, j, k), bufferVy_MinusX(j, k), vy(1, j, k));

                vs = Grid::MeanPhysicalProperty4(SwaveSpeed(nx - 1, j, k), SwaveSpeed(nx - 1, j + 1, k), SwaveSpeed(nx - 2, j, k), SwaveSpeed(nx - 2, j + 1, k));
                vDtDx = vs * dt / dx;
                vy(nx - 1, j, k) = Grid::UseNonreflectingCondition(vDtDx, spatialWeightNonReflective, temporalWeightNonReflective, vy(nx - 1, j, k), bufferVy_PlusX(j, k), vy(nx - 2, j, k));
            }
        }
#pragma omp for collapse(2)
        for (int j = 1; j < ny - 1; ++j)
        {
            for (int k = 0; k < nz - 1; ++k)
            {
                double vs = Grid::MeanPhysicalProperty4(SwaveSpeed(0, j, k), SwaveSpeed(1, j, k), SwaveSpeed(0, j, k + 1), SwaveSpeed(1, j, k + 1));
                double vDtDx = vs * dt / dx;
                vz(0, j, k) = Grid::UseNonreflectingCondition(vDtDx, spatialWeightNonReflective, temporalWeightNonReflective, vz(0, j, k), bufferVz_MinusX(j, k), vz(1, j, k));

                vs = Grid::MeanPhysicalProperty4(SwaveSpeed(nx - 1, j, k), SwaveSpeed(nx - 2, j, k), SwaveSpeed(nx - 1, j, k + 1), SwaveSpeed(nx - 2, j, k + 1));
                vDtDx = vs * dt / dx;
                vz(nx - 1, j, k) = Grid::UseNonreflectingCondition(vDtDx, spatialWeightNonReflective, temporalWeightNonReflective, vz(nx - 1, j, k), bufferVz_PlusX(j, k), vz(nx - 2, j, k));
            }
        }

        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        //
        // 2.Vx,VzをY方向に拡大
        //
        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#pragma omp for collapse(2)
        for (int k = 0; k < nz - 1; ++k)
        {
            for (int i = 0; i < nx - 1; ++i)
            {
                double vs = Grid::MeanPhysicalProperty4(SwaveSpeed(i, 0, k), SwaveSpeed(i, 1, k), SwaveSpeed(i + 1, 0, k), SwaveSpeed(i + 1, 1, k));
                double vDtDy = vs * dt / dy;
                vx(i, 0, k) = Grid::UseNonreflectingCondition(vDtDy, spatialWeightNonReflective, temporalWeightNonReflective, vx(i, 1, k), bufferVx_MinusY(k, i), vx(i, 1, k));

                vs = Grid::MeanPhysicalProperty4(SwaveSpeed(i, ny - 1, k), SwaveSpeed(i, ny - 2, k), SwaveSpeed(i + 1, ny - 1, k), SwaveSpeed(i + 1, ny - 2, k));
                vDtDy = vs * dt / dy;
                vx(i, ny - 1, k) = Grid::UseNonreflectingCondition(vDtDy, spatialWeightNonReflective, temporalWeightNonReflective, vx(i, ny - 1, k), bufferVx_PlusY(k, i), vx(i, ny - 2, k));
            }
        }
#pragma omp for collapse(2)
        for (int k = 0; k < nz - 1; ++k)
        {
            for (int i = 1; i < nx - 1; ++i)
            {
                double vs = Grid::MeanPhysicalProperty4(SwaveSpeed(i, 0, k), SwaveSpeed(i, 0, k + 1), SwaveSpeed(i, 1, k), SwaveSpeed(i, 1, k + 1));
                double vDtDy = vs * dt / dy;
                vz(i, 0, k) = Grid::UseNonreflectingCondition(vDtDy, spatialWeightNonReflective, temporalWeightNonReflective, vz(i, 1, k), bufferVz_MinusY(k, i), vz(i, 1, k));

                vs = Grid::MeanPhysicalProperty4(SwaveSpeed(i, ny - 1, k), SwaveSpeed(i, ny - 1, k + 1), SwaveSpeed(i, ny - 2, k), SwaveSpeed(i, ny - 2, k + 1));
                vDtDy = vs * dt / dy;
                vz(i, ny - 1, k) = Grid::UseNonreflectingCondition(vDtDy, spatialWeightNonReflective, temporalWeightNonReflective, vz(i, ny - 1, k), bufferVz_PlusY(k, i), vz(i, ny - 2, k));
            }
        }

        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        //
        // 3.VzをXY平面の四隅へ拡大
        //
        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#pragma omp for
        for (int k = 0; k < nz - 1; ++k)
        {
            double vsyz = Grid::MeanPhysicalProperty4(SwaveSpeed(0, 0, k), SwaveSpeed(0, 0, k + 1), SwaveSpeed(0, 1, k), SwaveSpeed(0, 1, k + 1));
            double vszx = Grid::MeanPhysicalProperty4(SwaveSpeed(0, 0, k), SwaveSpeed(1, 0, k), SwaveSpeed(0, 0, k + 1), SwaveSpeed(1, 0, k + 1));
            double vDtDy = vsyz * dt / dy;
            double vDtDx = vszx * dt / dx;
            vz(0, 0, k) = Grid::UseEdgeNonreflectingCondition(vDtDx, vDtDy, spatialWeightNonReflective, temporalWeightNonReflective,
                                                              vz(0, 0, k), bufferVz_MinusY(k, 0), bufferVz_MinusX(0, k), vz(0, 1, k), vz(1, 0, k));

            vsyz = Grid::MeanPhysicalProperty4(SwaveSpeed(0, ny - 1, k), SwaveSpeed(0, ny - 1, k + 1), SwaveSpeed(0, ny - 2, k), SwaveSpeed(0, ny - 2, k + 1));
            vszx = Grid::MeanPhysicalProperty4(SwaveSpeed(0, ny - 1, k), SwaveSpeed(1, ny - 1, k), SwaveSpeed(0, ny - 1, k + 1), SwaveSpeed(1, ny - 1, k + 1));
            vDtDy = vsyz * dt / dy;
            vDtDx = vszx * dt / dx;
            vz(0, ny - 1, k) = Grid::UseEdgeNonreflectingCondition(vDtDx, vDtDy, spatialWeightNonReflective, temporalWeightNonReflective,
                                                                   vz(0, ny - 1, k), bufferVz_PlusY(k, 0), bufferVz_MinusX(ny - 1, k), vz(0, ny - 2, k), vz(1, ny - 1, k));

            vsyz = Grid::MeanPhysicalProperty4(SwaveSpeed(nx - 1, 0, k), SwaveSpeed(nx - 1, 0, k + 1), SwaveSpeed(nx - 1, 1, k), SwaveSpeed(nx - 1, 1, k + 1));
            vszx = Grid::MeanPhysicalProperty4(SwaveSpeed(nx - 1, 0, k), SwaveSpeed(nx - 2, 0, k), SwaveSpeed(nx - 1, 0, k + 1), SwaveSpeed(nx - 2, 0, k + 1));
            vDtDy = vsyz * dt / dy;
            vDtDx = vszx * dt / dx;
            vz(nx - 1, 0, k) = Grid::UseEdgeNonreflectingCondition(vDtDx, vDtDy, spatialWeightNonReflective, temporalWeightNonReflective,
                                                                   vz(nx - 1, 0, k), bufferVz_MinusY(k, nx - 1), bufferVz_PlusX(0, k), vz(nx - 1, 1, k), vz(nx - 2, 0, k));

            vsyz = Grid::MeanPhysicalProperty4(SwaveSpeed(nx - 1, ny - 1, k), SwaveSpeed(nx - 1, ny - 1, k + 1), SwaveSpeed(nx - 1, ny - 2, k), SwaveSpeed(nx - 1, ny - 2, k + 1));
            vszx = Grid::MeanPhysicalProperty4(SwaveSpeed(nx - 1, ny - 1, k), SwaveSpeed(nx - 2, ny - 1, k), SwaveSpeed(nx - 1, ny - 1, k + 1), SwaveSpeed(nx - 2, ny - 1, k + 1));
            vDtDy = vsyz * dt / dy;
            vDtDx = vszx * dt / dx;
            vz(nx - 1, ny - 1, k) = Grid::UseEdgeNonreflectingCondition(vDtDx, vDtDy, spatialWeightNonReflective, temporalWeightNonReflective,
                                                                        vz(nx - 1, ny - 1, k), bufferVz_PlusY(k, nx - 1), bufferVz_PlusX(ny - 1, k), vz(nx - 1, ny - 2, k), vz(nx - 2, ny - 1, k));
        }

        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        //
        // 4.Vx,VyをZ方向に拡大
        //
        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#pragma omp for collapse(2)
        for (int i = 0; i < nx - 1; ++i)
        {
            for (int j = 0; j < ny; ++j)
            {
                double vs = Grid::MeanPhysicalProperty4(SwaveSpeed(i, j, nz - 1), SwaveSpeed(i + 1, j, nz - 1), SwaveSpeed(i, j, nz - 2), SwaveSpeed(i + 1, j, nz - 2));
                double vDtDz = vs * dt / dz;
                vx(i, j, nz - 1) = Grid::UseNonreflectingCondition(vDtDz, spatialWeightNonReflective, temporalWeightNonReflective, vx(i, j, nz - 1), bufferVx_PlusZ(i, j), vx(i, j, nz - 2));
            }
        }
#pragma omp for collapse(2)
        for (int i = 0; i < nx; ++i)
        {
            for (int j = 0; j < ny - 1; ++j)
            {
                double vs = Grid::MeanPhysicalProperty4(SwaveSpeed(i, j, nz - 1), SwaveSpeed(i, j, nz - 2), SwaveSpeed(i, j + 1, nz - 1), SwaveSpeed(i, j + 1, nz - 2));
                double vDtDz = vs * dt / dz;
                vy(i, j, nz - 1) = Grid::UseNonreflectingCondition(vDtDz, spatialWeightNonReflective, temporalWeightNonReflective, vy(i, j, nz - 1), bufferVy_PlusZ(i, j), vy(i, j, nz - 2));
            }
        }

        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        //
        // 5.VxをX方向に拡大
        //
        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#pragma omp for collapse(2)
        for (int j = 0; j < ny; ++j)
        {
            for (int k = 0; k < nz; ++k)
            {
                double vp = PwaveSpeed(0, j, k);
                double vDtDx = vp * dt / dx;
                vx(-1, j, k) = Grid::UseNonreflectingCondition(vDtDx, spatialWeightNonReflective, temporalWeightNonReflective, vx(-1, j, k), bufferVx_MinusX(j, k), vx(0, j, k));

                vp = PwaveSpeed(nx - 1, j, k);
                vDtDx = vp * dt / dx;
                vx(nx - 1, j, k) = Grid::UseNonreflectingCondition(vDtDx, spatialWeightNonReflective, temporalWeightNonReflective, vx(nx - 1, j, k), bufferVx_PlusX(j, k), vx(nx - 2, j, k));
            }
        }

        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        //
        // 6.VyをY方向に拡大
        //
        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#pragma omp for collapse(2)
        for (int k = 0; k < nz; ++k)
        {
            for (int i = 0; i < nx; ++i)
            {
                double vp = PwaveSpeed(i, 0, k);
                double vDtDy = vp * dt / dy;
                vy(i, -1, k) = Grid::UseNonreflectingCondition(vDtDy, spatialWeightNonReflective, temporalWeightNonReflective, vy(i, -1, k), bufferVy_MinusY(k, i), vy(i, 0, k));

                vp = PwaveSpeed(i, ny - 1, k);
                vDtDy = vp * dt / dy;
                vy(i, ny - 1, k) = Grid::UseNonreflectingCondition(vDtDy, spatialWeightNonReflective, temporalWeightNonReflective, vy(i, ny - 1, k), bufferVy_PlusY(k, i), vy(i, ny - 2, k));
            }
        }

        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        //
        // 7.VzをZ方向に拡大
        //
        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#pragma omp for collapse(2)
        for (int i = 0; i < nx; ++i)
        {
            for (int j = 0; j < ny; ++j)
            {
                double vp = PwaveSpeed(i, j, nz - 1);
                double vDtDz = vp * dt / dz;
                vz(i, j, nz - 1) = Grid::UseNonreflectingCondition(vDtDz, spatialWeightNonReflective, temporalWeightNonReflective, vz(i, j, nz - 1), bufferVz_PlusZ(i, j), vz(i, j, nz - 2));
            }
        }

        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        // 地表境界条件
        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#pragma omp for collapse(2)
        for (int i = 0; i < nx; ++i)
        {
            for (int j = 0; j < ny; ++j)
            {
                double lame1 = lambda(i, j, 0);
                double lame2 = mu(i, j, 0);
                double dxVx = DxVx(i, j, 0, useForthDerivative);
                double dyVy = DyVy(i, j, 0, useForthDerivative);
                double dzVz = -(dxVx + dyVy) * lame1 / (lame1 + 2.0 * lame2);

                vz(i, j, -1) = vz(i, j, 0) - dzVz * dz;
            }
        }
    }
}

void Grid::applyVelocityAbsorbingBoundary(double dt, int absorbingGridNumber, double absorbingDampingFactor)
{
    // Cerjan et al. 1985の吸収境界条件を適用

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // 吸収境界条件を適用
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // std::cout << "Applying absorbing boundary conditions" << std::endl;
    // std::cout << "  absorbingGridNumber = " << absorbingGridNumber << std::endl;
    // std::cout << "  absorbingDampingFactor = " << absorbingDampingFactor << std::endl;
    // std::cout << "" << std::endl;

    // 吸収境界条件の係数を計算
    Array1D<double> absorbingFactor;
    absorbingFactor.resize(0, absorbingGridNumber - 1);
    for (int i = 0; i < absorbingGridNumber; ++i)
    {
        absorbingFactor(i) = exp(-absorbingDampingFactor * absorbingDampingFactor * ((absorbingGridNumber - i) * (absorbingGridNumber - i)));
    }

#pragma omp parallel
    {

#pragma omp for collapse(3)
        for (int j = 0; j < ny; ++j)
        {
            for (int k = 0; k < nz; ++k)
            {
                for (int i = 0; i < absorbingFactor.size(); ++i)
                {
                    vx(i - 1, j, k) *= absorbingFactor(i);
                    vx(nx - 1 - i, j, k) *= absorbingFactor(i);
                }
            }
        }
#pragma omp for collapse(3)
        for (int j = -1; j < ny; ++j)
        {
            for (int k = 0; k < nz; ++k)
            {
                for (int i = 0; i < absorbingFactor.size(); ++i)
                {
                    vy(i, j, k) *= absorbingFactor(i);
                    vy(nx - 1 - i, j, k) *= absorbingFactor(i);
                }
            }
        }
#pragma omp for collapse(3)
        for (int j = 0; j < ny; ++j)
        {
            for (int k = -1; k < nz; ++k)
            {
                for (int i = 0; i < absorbingFactor.size(); ++i)
                {
                    vz(i, j, k) *= absorbingFactor(i);
                    vz(nx - 1 - i, j, k) *= absorbingFactor(i);
                }
            }
        }
#pragma omp for collapse(3)
        for (int k = 0; k < nz; ++k)
        {
            for (int i = -1; i < nx; ++i)
            {
                for (int j = 0; j < absorbingFactor.size(); ++j)
                {
                    vx(i, j, k) *= absorbingFactor(j);
                    vx(i, ny - 1 - j, k) *= absorbingFactor(j);
                }
            }
        }
#pragma omp for collapse(3)
        for (int k = 0; k < nz; ++k)
        {
            for (int i = 0; i < nx; ++i)
            {
                for (int j = 0; j < absorbingFactor.size(); ++j)
                {
                    vy(i, j - 1, k) *= absorbingFactor(j);
                    vy(i, ny - 1 - j, k) *= absorbingFactor(j);
                }
            }
        }
#pragma omp for collapse(3)
        for (int k = -1; k < nz; ++k)
        {
            for (int i = 0; i < nx; ++i)
            {
                for (int j = 0; j < absorbingFactor.size(); ++j)
                {
                    vz(i, j, k) *= absorbingFactor(j);
                    vz(i, ny - 1 - j, k) *= absorbingFactor(j);
                }
            }
        }
#pragma omp for collapse(3)
        for (int i = -1; i < nx; ++i)
        {
            for (int j = 0; j < ny; ++j)
            {
                for (int k = 0; k < absorbingFactor.size(); ++k)
                {
                    vx(i, j, nz - 1 - k) *= absorbingFactor(k);
                }
            }
        }
#pragma omp for collapse(3)
        for (int i = 0; i < nx; ++i)
        {
            for (int j = -1; j < ny; ++j)
            {
                for (int k = 0; k < absorbingFactor.size(); ++k)
                {
                    vy(i, j, nz - 1 - k) *= absorbingFactor(k);
                }
            }
        }
#pragma omp for collapse(3)
        for (int i = 0; i < nx; ++i)
        {
            for (int j = 0; j < ny; ++j)
            {
                for (int k = 0; k < absorbingFactor.size(); ++k)
                {
                    vz(i, j, nz - 1 - k) *= absorbingFactor(k);
                }
            }
        }
    }
}

void Grid::updateStress(double dt, double f0, bool useForthDerivative)
{
    // 応力更新式に基づいて、各グリッド点の応力成分を計算
    // ... (Graves, 1996およびPitarka, 1999の応力更新式を参照)

    // std::cout << "updateStress" << std::endl;
    // std::cout << "" << std::endl;

#pragma omp parallel
    {
        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        // Txx, Tyy, Tzzの更新
        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#pragma omp for collapse(3)
        for (int i = 0; i < nx; ++i)
        {
            for (int j = 0; j < ny; ++j)
            {
                for (int k = 0; k < nz; ++k)
                {
                    double lame1 = lambda(i, j, k);
                    double lame2 = mu(i, j, k);
                    double Q0val = Q0(i, j, k);
                    double damping = Grid::CalculateDampingFactor(f0, dt, Q0val);

                    double dxVx = DxVx(i, j, k, useForthDerivative);
                    double dyVy = DyVy(i, j, k, useForthDerivative);
                    double dzVz = DzVz(i, j, k, useForthDerivative);

                    txx(i, j, k) = damping * (txx(i, j, k) + (lame1 * (dxVx + dyVy + dzVz) + 2.0 * lame2 * dxVx) * dt);
                    tyy(i, j, k) = damping * (tyy(i, j, k) + (lame1 * (dxVx + dyVy + dzVz) + 2.0 * lame2 * dyVy) * dt);
                    if (k == 0)
                    {
                        tzz(i, j, k) = 0.0; // z=0 (境界条件)
                    }
                    else
                    {
                        tzz(i, j, k) = damping * (tzz(i, j, k) + (lame1 * (dxVx + dyVy + dzVz) + 2.0 * lame2 * dzVz) * dt);
                    }
                }
            }
        }
        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        // Txy, Txz, Tyzの更新
        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#pragma omp for collapse(3)
        for (int i = 0; i < nx; ++i)
        {
            for (int j = 0; j < ny - 1; ++j)
            {
                for (int k = 0; k < nz - 1; ++k)
                {
                    double lame2 = Grid::MeanPhysicalProperty4(mu(i, j, k), mu(i, j, k + 1), mu(i, j + 1, k), mu(i, j + 1, k + 1));
                    double Q0val = Grid::MeanPhysicalProperty4(Q0(i, j, k), Q0(i, j, k + 1), Q0(i, j + 1, k), Q0(i, j + 1, k + 1));
                    double damping = Grid::CalculateDampingFactor(f0, dt, Q0val);
                    double dzVy = DzVy(i, j, k, useForthDerivative);
                    double dyVz = DyVz(i, j, k, useForthDerivative);
                    tyz(i, j, k) = (damping * (tyz(i, j, k) + lame2 * (dzVy + dyVz) * dt));
                }
            }
        }
#pragma omp for collapse(3)
        for (int i = 0; i < nx - 1; ++i)
        {
            for (int j = 0; j < ny; ++j)
            {
                for (int k = 0; k < nz - 1; ++k)
                {
                    double lame2 = Grid::MeanPhysicalProperty4(mu(i, j, k), mu(i + 1, j, k), mu(i, j, k + 1), mu(i + 1, j, k + 1));
                    double Q0val = Grid::MeanPhysicalProperty4(Q0(i, j, k), Q0(i + 1, j, k), Q0(i, j, k + 1), Q0(i + 1, j, k + 1));
                    double damping = Grid::CalculateDampingFactor(f0, dt, Q0val);
                    double dxVz = DxVz(i, j, k, useForthDerivative);
                    double dzVx = DzVx(i, j, k, useForthDerivative);
                    tzx(i, j, k) = (damping * (tzx(i, j, k) + lame2 * (dxVz + dzVx) * dt));
                }
            }
        }
#pragma omp for collapse(3)
        for (int i = 0; i < nx - 1; ++i)
        {
            for (int j = 0; j < ny - 1; ++j)
            {
                for (int k = 0; k < nz; ++k)
                {
                    double lame2 = Grid::MeanPhysicalProperty4(mu(i, j, k), mu(i, j + 1, k), mu(i + 1, j, k), mu(i + 1, j + 1, k));
                    double Q0val = Grid::MeanPhysicalProperty4(Q0(i, j, k), Q0(i, j + 1, k), Q0(i + 1, j, k), Q0(i + 1, j + 1, k));
                    double damping = Grid::CalculateDampingFactor(f0, dt, Q0val);
                    double dyVx = DyVx(i, j, k, useForthDerivative);
                    double dxVy = DxVy(i, j, k, useForthDerivative);
                    txy(i, j, k) = (damping * (txy(i, j, k) + lame2 * (dyVx + dxVy) * dt));
                }
            }
        }
    }
}

void Grid::applyStressAbsorbingBoundary(double dt, int absorbingGridNumber, double absorbingDampingFactor)
{
    // Cerjan et al. 1985の吸収境界条件を適用

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // 吸収境界条件を適用
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // std::cout << "Applying absorbing boundary conditions" << std::endl;
    // std::cout << "  absorbingGridNumber = " << absorbingGridNumber << std::endl;
    // std::cout << "  absorbingDampingFactor = " << absorbingDampingFactor << std::endl;
    // std::cout << "" << std::endl;

    // 吸収境界条件の係数を計算
    Array1D<double> absorbingFactor;
    absorbingFactor.resize(0, absorbingGridNumber - 1);
    for (int i = 0; i < absorbingGridNumber; ++i)
    {
        absorbingFactor(i) = exp(-absorbingDampingFactor * absorbingDampingFactor * (absorbingGridNumber - i) * (absorbingGridNumber - i));
    }

#pragma omp parallel
    {
#pragma omp for collapse(3)
        for (int j = 0; j < ny; ++j)
        {
            for (int k = 0; k < nz; ++k)
            {
                for (int i = 0; i < absorbingFactor.size(); ++i)
                {
                    txx(i, j, k) *= absorbingFactor(i);
                    tyy(i, j, k) *= absorbingFactor(i);
                    tzz(i, j, k) *= absorbingFactor(i);
                    txx(nx - 1 - i, j, k) *= absorbingFactor(i);
                    tyy(nx - 1 - i, j, k) *= absorbingFactor(i);
                    tzz(nx - 1 - i, j, k) *= absorbingFactor(i);
                }
            }
        }
#pragma omp for collapse(3)
        for (int j = 0; j < ny - 1; ++j)
        {
            for (int k = 0; k < nz - 1; ++k)
            {
                for (int i = 0; i < absorbingFactor.size(); ++i)
                {
                    tyz(i, j, k) *= absorbingFactor(i);
                    tyz(nx - 1 - i, j, k) *= absorbingFactor(i);
                }
            }
        }
#pragma omp for collapse(3)
        for (int j = 0; j < ny; ++j)
        {
            for (int k = 0; k < nz - 1; ++k)
            {
                for (int i = 0; i < absorbingFactor.size(); ++i)
                {
                    tzx(i, j, k) *= absorbingFactor(i);
                    tzx(nx - 2 - i, j, k) *= absorbingFactor(i);
                }
            }
        }
#pragma omp for collapse(3)
        for (int j = 0; j < ny - 1; ++j)
        {
            for (int k = 0; k < nz; ++k)
            {
                for (int i = 0; i < absorbingFactor.size(); ++i)
                {
                    txy(i, j, k) *= absorbingFactor(i);
                    txy(nx - 2 - i, j, k) *= absorbingFactor(i);
                }
            }
        }
#pragma omp for collapse(3)
        for (int k = 0; k < nz; ++k)
        {
            for (int i = 0; i < nx; ++i)
            {
                for (int j = 0; j < absorbingFactor.size(); ++j)
                {
                    txx(i, j, k) *= absorbingFactor(j);
                    tyy(i, j, k) *= absorbingFactor(j);
                    tzz(i, j, k) *= absorbingFactor(j);
                    txx(i, ny - 1 - j, k) *= absorbingFactor(j);
                    tyy(i, ny - 1 - j, k) *= absorbingFactor(j);
                    tzz(i, ny - 1 - j, k) *= absorbingFactor(j);
                }
            }
        }
#pragma omp for collapse(3)
        for (int k = 0; k < nz - 1; ++k)
        {
            for (int i = 0; i < nx; ++i)
            {
                for (int j = 0; j < absorbingFactor.size(); ++j)
                {
                    tyz(i, j, k) *= absorbingFactor(j);
                    tyz(i, ny - 2 - j, k) *= absorbingFactor(j);
                }
            }
        }
#pragma omp for collapse(3)
        for (int k = 0; k < nz - 1; ++k)
        {
            for (int i = 0; i < nx - 1; ++i)
            {
                for (int j = 0; j < absorbingFactor.size(); ++j)
                {
                    tzx(i, j, k) *= absorbingFactor(j);
                    tzx(i, ny - 1 - j, k) *= absorbingFactor(j);
                }
            }
        }
#pragma omp for collapse(3)
        for (int k = 0; k < nz; ++k)
        {
            for (int i = 0; i < nx - 1; ++i)
            {
                for (int j = 0; j < absorbingFactor.size(); ++j)
                {
                    txy(i, j, k) *= absorbingFactor(j);
                    txy(i, ny - 2 - j, k) *= absorbingFactor(j);
                }
            }
        }
#pragma omp for collapse(3)
        for (int i = 0; i < nx; ++i)
        {
            for (int j = 0; j < ny; ++j)
            {
                for (int k = 0; k < absorbingFactor.size(); ++k)
                {
                    txx(i, j, nz - 1 - k) *= absorbingFactor(k);
                    tyy(i, j, nz - 1 - k) *= absorbingFactor(k);
                    tzz(i, j, nz - 1 - k) *= absorbingFactor(k);
                }
            }
        }
#pragma omp for collapse(3)
        for (int i = 0; i < nx; ++i)
        {
            for (int j = 0; j < ny - 1; ++j)
            {
                for (int k = 0; k < absorbingFactor.size(); ++k)
                {
                    tyz(i, j, nz - 2 - k) *= absorbingFactor(k);
                }
            }
        }
#pragma omp for collapse(3)
        for (int i = 0; i < nx - 1; ++i)
        {
            for (int j = 0; j < ny; ++j)
            {
                for (int k = 0; k < absorbingFactor.size(); ++k)
                {
                    tzx(i, j, nz - 2 - k) *= absorbingFactor(k);
                }
            }
        }
#pragma omp for collapse(3)
        for (int i = 0; i < nx - 1; ++i)
        {
            for (int j = 0; j < ny - 1; ++j)
            {
                for (int k = 0; k < absorbingFactor.size(); ++k)
                {
                    txy(i, j, nz - 1 - k) *= absorbingFactor(k);
                }
            }
        }
    }
}

void Grid::applySource(double dt, int i, int j, int k, double momentTensor[3][3])
{
    // momentTensor[3][3] : 震源時間関数　モーメントテンソルを1階時間微分したもの
    // ... (応力成分に震源励源項を加える)

    // std::cout << "applySource" << std::endl;
    // std::cout << "" << std::endl;

    // index で指定する　iiはi,j,k ijは隣接格子に分配

    double V = dx * dy * dz; // 均等格子の場合 dx * dy * dz

    double mxx = dt / V * momentTensor[0][0];
    double myy = dt / V * momentTensor[1][1];
    double mzz = dt / V * momentTensor[2][2];
    double myz = dt / V / 4.0 * momentTensor[1][2];
    double mzx = dt / V / 4.0 * momentTensor[2][0];
    double mxy = dt / V / 4.0 * momentTensor[0][1];

    txx(i, j, k) -= mxx;
    tyy(i, j, k) -= myy;
    tzz(i, j, k) -= mzz;

    tyz(i, j, k) -= myz; // myz
    tyz(i, j, k - 1) -= myz;
    tyz(i, j - 1, k) -= myz;
    tyz(i, j - 1, k - 1) -= myz;

    tzx(i, j, k) -= mzx;
    tzx(i, j, k - 1) -= mzx;
    tzx(i - 1, j, k) -= mzx;
    tzx(i - 1, j, k - 1) -= mzx;

    txy(i, j, k) -= mxy;
    txy(i, j - 1, k) -= mxy;
    txy(i - 1, j, k) -= mxy;
    txy(i - 1, j - 1, k) -= mxy;

    // std::cout << "applySource" << std::endl;
    // std::cout << "mxx = " << mxx << std::endl;
    // std::cout << "myy = " << myy << std::endl;
    // std::cout << "mzz = " << mzz << std::endl;
    // std::cout << "myz = " << myz << std::endl;
    // std::cout << "mzx = " << mzx << std::endl;
    // std::cout << "mxy = " << mxy << std::endl;
    // std::cout << "" << std::endl;
}

void Grid::getVelocityAtReceiver(int i, int j, int k, double &ux, double &uy, double &uz) const
{
    // 受信点での速度成分を取得
    ux = vx(i, j, k);
    uy = vy(i, j, k);
    uz = vz(i, j, k);
}

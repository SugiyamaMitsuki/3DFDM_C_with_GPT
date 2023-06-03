#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include <iomanip>
#include <omp.h> // OpenMP

#include "array3d.h"
#include "Simulation.h"
#include "CoordinateManager.h"
#include "Config.h"
#include "SourceInfo.h"
#include "ReceiverInfo.h"
#include "Grid.h"
#include "Source.h"
#include "SlipVelocityInfo.h"
#include "Receiver.h"

Simulation::Simulation(const std::string &configFile)
    : grid()
{
    readConfig(configFile);

    omp_set_num_threads(num_threads);
    std::cout << "Simulation parameters:" << std::endl
              << "  nx: " << nx << std::endl
              << "  ny: " << ny << std::endl
              << "  nz: " << nz << std::endl
              << "  dx: " << dx << std::endl
              << "  dy: " << dy << std::endl
              << "  dz: " << dz << std::endl
              << "  dt: " << dt << std::endl
              << "  totalTime: " << totalTime << std::endl
              << "  f0: " << f0 << std::endl
              << "  spatialWeightNonReflective: " << spatialWeightNonReflective << std::endl
              << "  temporalWeightNonReflective: " << temporalWeightNonReflective << std::endl
              << "  absorbingGridNumber: " << absorbingGridNumber << std::endl
              << "  absorbingDampingFactor: " << absorbingDampingFactor << std::endl
              << std::endl;
}

void Simulation::readConfig(const std::string &configFile)
{
    // 設定ファイルの読み込み
    Config config(configFile);

    // 入力パラメータの設定
    nx = config.get<int>("nx");
    ny = config.get<int>("ny");
    nz = config.get<int>("nz");
    dx = config.get<double>("dx");
    dy = config.get<double>("dy");
    dz = config.get<double>("dz");
    totalTime = config.get<double>("totalTime");
    dt = config.get<double>("dt");
    snapshotInterval = config.get<int>("snapshotInterval");
    num_threads = config.get<int>("num_threads");
    f0 = config.get<double>("f0");
    spatialWeightNonReflective = config.get<double>("spatialWeightNonReflective");
    temporalWeightNonReflective = config.get<double>("temporalWeightNonReflective");
    absorbingGridNumber = config.get<int>("absorbingGridNumber");
    absorbingDampingFactor = config.get<double>("absorbingDampingFactor");
    lat0 = config.get<double>("lat0");
    lon0 = config.get<double>("lon0");
    depth0 = config.get<double>("depth0") * 1000.0;
}

void Simulation::initializeModel()
{
    grid.initializeModel(nx, ny, nz, dx, dy, dz, lat0, lon0, depth0);
}
void Simulation::setPhysicalProperty(bool isTest, std::string inputPhysicalProperty, std::string inputLayerProperty)
{
    grid.setPhysicalProperty(isTest, inputPhysicalProperty, inputLayerProperty);
}
void Simulation::addSource(const Source &source)
{
    sources.push_back(source);
}

void Simulation::setSource()
{
    std::vector<SourceInfo> sourceInfo = readSourcesCSV("./input/sources.csv");
    for (const auto &info : sourceInfo)
    {
        // id = 0 はダミー
        if (info.id == 0)
        {
            continue;
        }
        std::tuple<int, int, int> gridIndex = CoordinateManager::LatLonDepthToIndex(info.lat, info.lon, info.depth, nx, ny, nz, dx, dy, dz, lat0, lon0, depth0);
        int xgrid = std::get<0>(gridIndex);
        int ygrid = std::get<1>(gridIndex);
        int zgrid = std::get<2>(gridIndex);
        Source source1(xgrid, ygrid, zgrid, info.M0, info.strike, info.dip, info.rake);
        // 震源位置はgridのindexで指定してるけど、実際は座標で指定したい　近隣のgridに重みをつけて分配したい
        Simulation::addSource(source1);

        // debug
        //  std::cout << "Source" << info.id << ":" << std::endl;
        //  std::cout << "  xgrid: " << xgrid << std::endl;
        //  std::cout << "  ygrid: " << ygrid << std::endl;
        //  std::cout << "  zgrid: " << zgrid << std::endl;
        //  std::cout << "  M0: " << info.M0 << std::endl;
        //  std::cout << "  strike: " << info.strike << std::endl;
        //  std::cout << "  dip: " << info.dip << std::endl;
        //  std::cout << "  rake: " << info.rake << std::endl;
    }

    // slipVelocityの読み込み
    SlipVelocityReader slipVelocityReader("./input/slipVelocity.csv");
    std::cout << "slipVelocityReader" << std::endl;
    const auto &slipVelocitys = slipVelocityReader.getData();

    // setしたsourceに、slipVelocityを割り当てる
    int count = -1;
    for (auto &src : sources)
    {
        // debug
        // std::cout << "Source" << count + 1 << ":" << std::endl;
        count++;
        const auto &slipVelocity = slipVelocitys[count];
        Array1D<double> inputTime(slipVelocity.time);
        Array1D<double> inputSlipVelocity(slipVelocity.velocity);

        src.setSlipVelocity(inputTime, inputSlipVelocity, dt, totalTime);
    }

    // Moment-rate time history 中村・宮武、地震-第2輯vol53-No1
    // double Vm = 1.7;      // maximum slip velocity (m/s)
    // double td = 0.05;     // max slip velocity time (s) fmax=1/(pi*td)
    // double tb = 0.08;     // transition time (s)
    // double tr = 0.7;      // rise time (s)
    // double ts = 1.5 * tr; // slip end time (s)
    // src.makeSlipVelocity(Vm, td, tb, tr, ts, dt, totalTime);
}

void Simulation::addReceiver(const Receiver &receiver)
{
    receivers.push_back(receiver);
}
void Simulation::setReceiver()
{
    std::vector<ReceiverInfo> receiverInfo = readReceiversCSV("./input/receivers.csv");
    for (const auto &rec : receiverInfo)
    {
        // id = 0 はダミー
        if (rec.id == 0)
        {
            continue;
        }
        std::tuple<int, int, int> gridIndex = CoordinateManager::LatLonDepthToIndex(rec.lat, rec.lon, rec.depth, nx, ny, nz, dx, dy, dz, lat0, lon0, depth0);
        int xgrid = std::get<0>(gridIndex);
        int ygrid = std::get<1>(gridIndex);
        int zgrid = std::get<2>(gridIndex);
        Receiver receiver1(rec.id, rec.lat, rec.lon, rec.depth, rec.name, xgrid, ygrid, zgrid);
        // 震源位置はgridのindexで指定してるけど、実際は座標で指定したい　近隣のgridに重みをつけて分配したい
        addReceiver(receiver1);

        // // debug
        // std::cout << "Receiver: " << "id" << rec.id << ", " << "xgrid" << xgrid << ", " << "ygrid" << ygrid << ", " << "zgrid" << zgrid << std::endl;
        // std::cout << "Receiver: " << "id" << rec.id << ", " << "lat" << rec.lat << ", " << "lon" << rec.lon << ", " << "depth" << rec.depth << std::endl;
        // std::cout << "" << std::endl;
    }
}
// Main loop
void Simulation::run()
{
    std::cout << "Simulation started." << std::endl;
    // シミュレーションの実行
    std::cout << "Simulation start." << std::endl;
    std::cout << "" << std::endl;

    std::cout << "Initialize model..." << std::endl;
    initializeModel();
    std::cout << "Initialize model done." << std::endl;
    std::cout << "" << std::endl;

    std::cout << "Set source..." << std::endl;
    setSource();
    std::cout << "Set source done." << std::endl;
    std::cout << "" << std::endl;

    std::cout << "Set receiver..." << std::endl;
    setReceiver();
    std::cout << "Set receiver done." << std::endl;
    std::cout << "" << std::endl;

    num_steps = totalTime / dt;

    bool useForthDerivative = true; // ハードコーディング　後で修正

    std::cout << "Set physical property..." << std::endl;
    bool isTest = false;
    std::string inputPhysicalProperty = "./GroundStructureData/D-V3.2-STRUCT_DEEP-PYS.csv"; // ハードコーディング　後で修正
    std::string inputLayerProperty = "./GroundStructureData/D-V3.2-STRUCT_DEEP-LYRE.csv";   // J-SHS 標高データ
    setPhysicalProperty(isTest, inputPhysicalProperty, inputLayerProperty);
    std::cout << "Set physical property done." << std::endl;
    std::cout << "" << std::endl;

    std::cout << "Set boundary condition..." << std::endl;
    std::string outputPhysicalProperty = "./output/Output_physicalProperty.csv";
    savePhysicalProperty(outputPhysicalProperty);
    std::cout << "Set boundary condition done." << std::endl;
    std::cout << "" << std::endl;

    std::cout << "Start simulation..." << std::endl;
    std::cout << "" << std::endl;
    for (int time_step = 0; time_step < num_steps; ++time_step)
    {
        std::cout << "Time step: " << time_step << " / " << num_steps << " (" << time_step * dt << " sec)" << std::endl;

        double VelocityTime = dt * (time_step - 1);
        double StressTime = dt * (time_step - 0.5);

        // 速度ベクトルを更新
        grid.updateVelocity(dt, useForthDerivative, f0, spatialWeightNonReflective, temporalWeightNonReflective);

        // モーメントでなく速度を指定する場合はここ 未実装
        if (false)
        {
            std::cout << "Not implemented yet." << std::endl;
        }

        // 吸収境界条件を適用
        grid.applyVelocityAbsorbingBoundary(dt, absorbingGridNumber, absorbingDampingFactor);

        // 応力テンソルを更新
        grid.updateStress(dt, f0, useForthDerivative);

        // 震源モーメントテンソルを更新
        applySourceAtTimeStep(dt, time_step); // 速度or応力テンソルによる指定を選択できるようにしておきたい 　　　　　　　　　　　　　　　　　震源時間関数を作成する　近隣グリッドに応力分配未実装

        // 吸収境界条件を適用
        grid.applyStressAbsorbingBoundary(dt, absorbingGridNumber, absorbingDampingFactor);

        // 受信点の波形を記録
        recordReceiversAtTimeStep(time_step);

        // スナップショットを記録
        if (time_step % snapshotInterval == 0)
        {
            std::cout << "Recording snapshot..." << std::endl;
            recordSnapshot(time_step);
            std::cout << "Recording snapshot done." << std::endl;
        }

        std::cout << "" << std::endl;
    }
    std::cout << "Simulation done." << std::endl;
    std::cout << "" << std::endl;

    std::cout << "Recording receiver data..." << std::endl;
    writeData("./output/Output_simulation.csv"); // Receiverのデータを書き出す
    std::cout << "Recording receiver data done." << std::endl;
    std::cout << "" << std::endl;

    std::cout << "Simulation finished." << std::endl;
}

void Simulation::savePhysicalProperty(std::string outputPhysicalProperty)
{
    std::cout << "Recording PhysicalProperty..." << std::endl;

    Array3D<double> &PwaveSpeedRef = grid.getPwaveSpeed();
    Array3D<double> &SwaveSpeedRef = grid.getSwaveSpeed();
    Array3D<double> &rhoRef = grid.getRho();
    Array3D<double> &Q0Ref = grid.getQ0();
    Array3D<double> &lambdaRef = grid.getLambda();
    Array3D<double> &muRef = grid.getMu();

    // ファイルをバイナリモードでオープン
    std::ofstream ofs(outputPhysicalProperty, std::ios::binary);
    if (!ofs)
    {
        std::cerr << "Error opening file for snapshot: " << outputPhysicalProperty << std::endl;
        return;
    }
    // C++のストリームとCのストリームの同期を無効化
    std::ios::sync_with_stdio(false);

    // 文字列バッファを用意
    std::ostringstream ss;

    // header
    ss << "x,y,z,lat,lon,depth,PwaveSpeed,SwaveSpeed,rho,Q0,lambda,mu" << std::endl;

#pragma omp parallel
    {
        std::ostringstream local_ss;

#pragma omp for schedule(dynamic)
        for (int i = 0; i < nx - 1; ++i)
        {
            for (int j = 0; j < ny - 1; ++j)
            {
                for (int k = 0; k < nz - 1; ++k)
                {
                    // 0~10層 debug用
                    if (k > 10)
                    {
                        continue;
                    }
                    std::tuple<double, double, double> latlondepth = CoordinateManager::IndexToLatLonDepth(i, j, k, nx, ny, nz, dx, dy, dz, lat0, lon0, depth0);
                    double lat = std::get<0>(latlondepth);
                    double lon = std::get<1>(latlondepth);
                    double depth = std::get<2>(latlondepth);
                    local_ss << i << "," << j << "," << k << "," << lat << "," << lon << "," << depth << "," << PwaveSpeedRef(i, j, k) << "," << SwaveSpeedRef(i, j, k) << "," << rhoRef(i, j, k) << "," << Q0Ref(i, j, k) << "," << lambdaRef(i, j, k) << "," << muRef(i, j, k) << std::endl;
                }
            }
        }
#pragma omp critical
        ss << local_ss.str();
    }
    // 文字列バッファの内容を一度にファイルに書き込む
    ofs << ss.str();
    ofs.close();
}

void Simulation::recordSnapshot(int time_step)
{
    Array3D<double> &vxRef = grid.getVx();
    Array3D<double> &vyRef = grid.getVy();
    Array3D<double> &vzRef = grid.getVz();

    std::string snapshotFilename = "./output/snapshot_";
    std::string timeStepString = std::to_string(time_step);
    std::string paddedTimeStep = std::string(5 - timeStepString.length(), '0') + timeStepString;
    std::string filename = snapshotFilename + paddedTimeStep + ".bin";

    std::vector<double> data_lat, data_lon, data_depth, data_vx, data_vy, data_vz;

#pragma omp parallel
    {
        std::vector<double> local_lat, local_lon, local_depth, local_vx, local_vy, local_vz;

#pragma omp for schedule(dynamic)
        for (int i = 0; i < nx - 1; ++i)
        {
            for (int j = 0; j < ny - 1; ++j)
            {
                for (int k = 0; k < nz - 1; ++k)
                {
                    if (k == 0) // test　表層だけ保存
                    {
                        std::tuple<double, double, double> latlondepth = CoordinateManager::IndexToLatLonDepth(i, j, k, nx, ny, nz, dx, dy, dz, lat0, lon0, depth0);
                        local_lat.push_back(std::get<0>(latlondepth));
                        local_lon.push_back(std::get<1>(latlondepth));
                        local_depth.push_back(std::get<2>(latlondepth));
                        local_vx.push_back(vxRef(i, j, k));
                        local_vy.push_back(vyRef(i, j, k));
                        local_vz.push_back(vzRef(i, j, k));
                    }
                }
            }
        }

#pragma omp critical
        {
            data_lat.insert(data_lat.end(), local_lat.begin(), local_lat.end());
            data_lon.insert(data_lon.end(), local_lon.begin(), local_lon.end());
            data_depth.insert(data_depth.end(), local_depth.begin(), local_depth.end());
            data_vx.insert(data_vx.end(), local_vx.begin(), local_vx.end());
            data_vy.insert(data_vy.end(), local_vy.begin(), local_vy.end());
            data_vz.insert(data_vz.end(), local_vz.begin(), local_vz.end());
        }
    }
    std::ofstream ofs(filename, std::ios::binary);
    if (!ofs)
    {
        std::cerr << "Error opening file for snapshot: " << filename << std::endl;
        return;
    }
    ofs.write(reinterpret_cast<char *>(data_lat.data()), data_lat.size() * sizeof(double));
    ofs.write(reinterpret_cast<char *>(data_lon.data()), data_lon.size() * sizeof(double));
    ofs.write(reinterpret_cast<char *>(data_depth.data()), data_depth.size() * sizeof(double));
    ofs.write(reinterpret_cast<char *>(data_vx.data()), data_vx.size() * sizeof(double));
    ofs.write(reinterpret_cast<char *>(data_vy.data()), data_vy.size() * sizeof(double));
    ofs.write(reinterpret_cast<char *>(data_vz.data()), data_vz.size() * sizeof(double));

    ofs.close();
}

void Simulation::applySourceAtTimeStep(double dt, int time_step)
{
    for (auto &src : sources)
    {
        double scale = src.getSlipVelocityAtTimeStep(time_step);

        double momentTensor[3][3];
        src.calculateMomentTensor(momentTensor, scale);

        // 震源位置を取得
        int i, j, k;
        src.getLocation(i, j, k);

        // std::cout << "    applySourceAtTimeStep Source: " << i << ", " << j << ", " << k << std::endl; // debug

        grid.applySource(dt, i, j, k, momentTensor);

        // 震源モーメントテンソルを表示 debug
        // std::cout << "        Source: " << i << ", " << j << ", " << k << std::endl;
        // std::cout << "        Moment tensor: " << std::endl;
        // std::cout << "            " << momentTensor[0][0] << ", " << momentTensor[0][1] << ", " << momentTensor[0][2] << std::endl;
        // std::cout << "            " << momentTensor[1][0] << ", " << momentTensor[1][1] << ", " << momentTensor[1][2] << std::endl;
        // std::cout << "            " << momentTensor[2][0] << ", " << momentTensor[2][1] << ", " << momentTensor[2][2] << std::endl;
        // std::cout << "" << std::endl;
    }
}

void Simulation::recordReceiversAtTimeStep(int time_step)
{
    // std::cout << "Recording receivers..." << std::endl;

    for (auto &receiver : receivers)
    {
        // 受信点位置を取得
        int x, y, z;
        receiver.getLocation(x, y, z);
        double ux, uy, uz;
        grid.getVelocityAtReceiver(x, y, z, ux, uy, uz);
        receiver.record(ux, uy, uz);
        // 受信波形データを表示 debug
        // std::cout << "        Receiver: " << x << ", " << y << ", " << z << std::endl;
        // std::cout << "        ux: " << ux << ", uy: " << uy << ", uz: " << uz << std::endl;
        // std::cout << "" << std::endl;
    }
}

void Simulation::writeData(const std::string &baseFilename)
{
    std::cout << "Writing data..." << std::endl;

    for (auto &receiver : receivers)
    {
        // 観測点名
        int id;
        double lat;
        double lon;
        double depth;
        std::string receiverName;
        receiver.getReceiverInfo(id, lat, lon, depth,receiverName);

        std::cout << "    Receiver: " << receiverName << std::endl;

        // baseFilename_receiverName_lon=?_lat=?_depth=?_ux,uy,uz.csv
        std::string filename = baseFilename + "_" + receiverName + "_lon=" + std::to_string(lon) + "_lat=" + std::to_string(lat) + "_depth=" + std::to_string(depth) + "_ux,uy,uz.csv";
        std::ofstream ofs(filename);
        if (!ofs)
        {
            std::cerr << "Cannot open file: " << filename << std::endl;
            return;
        }
        std::vector<double> uxData, uyData, uzData;
        receiver.getVelocityData(uxData, uyData, uzData);
        for (int i = 0; i < uxData.size(); ++i)
        {
            ofs << uxData[i] << "," << uyData[i] << "," << uzData[i] << std::endl;
            // std::cout << uxData[i] << "," << uyData[i] << "," << uzData[i] << std::endl; // debug
        }
        ofs.close();
    }
}

// void Simulation::recordSnapshot(int time_step)
// {
//     std::cout << "Recording snapshot..." << std::endl;

//     Array3D<double> &vxRef = grid.getVx();
//     Array3D<double> &vyRef = grid.getVy();
//     Array3D<double> &vzRef = grid.getVz();
//     // 出力ファイル名
//     std::string snapshotFilename = "./output/snapshot";
//     std::string filename = snapshotFilename + std::to_string(time_step) + ".csv";

//     // ファイルをバイナリモードでオープン
//     std::ofstream ofs(filename, std::ios::binary);
//     if (!ofs)
//     {
//         std::cerr << "Error opening file for snapshot: " << filename << std::endl;
//         return;
//     }

//     // C++のストリームとCのストリームの同期を無効化
//     std::ios::sync_with_stdio(false);

//     // 文字列バッファを用意
//     std::ostringstream ss;

//     // header
//     ss << "x,y,z,lat,lon,depth,vx,vy,vz" << std::endl;

//     #pragma omp parallel
//     {
//         std::ostringstream local_ss;

//         #pragma omp for schedule(dynamic)
//         for (int i = 0; i < nx - 1; ++i)
//         {
//             for (int j = 0; j < ny - 1; ++j)
//             {
//                 for (int k = 0; k < nz - 1; ++k)
//                 {
//                     if (k == 0)//test　表層だけ保存
//                     {
//                     std::tuple<double, double, double> latlondepth = CoordinateManager::IndexToLatLonDepth(i, j, k, nx, ny, nz, dx, dy, dz, lat0, lon0, depth0);
//                     double lat = std::get<0>(latlondepth);
//                     double lon = std::get<1>(latlondepth);
//                     double depth = std::get<2>(latlondepth);
//                     local_ss << i << "," << j << "," << k << "," << lat << "," << lon << "," << depth << "," << vxRef(i, j, k) << "," << vyRef(i, j, k) << "," << vzRef(i, j, k) << std::endl;
//                     }
//                 }
//             }
//         }

//         #pragma omp critical
//         ss << local_ss.str();
//     }
//     // 文字列バッファの内容を一度にファイルに書き込む
//     ofs << ss.str();
//     ofs.close();
// }
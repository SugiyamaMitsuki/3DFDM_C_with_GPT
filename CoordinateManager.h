#include <iostream>
#include <cmath>
constexpr double EARTH_RADIUS = 6371e3; // 地球の半径 (m)
// # 緯度経度の縮尺比調整
constexpr double Rx = 2.0 * M_PI * 6378137 / 1000.0;           // #赤道周囲 km　https://ja.wikipedia.org/wiki/GRS80
constexpr double Ry = 2.0 * M_PI * 6356752.314140356 / 1000.0; // #子午線周囲 km

class CoordinateManager
{
public:
    CoordinateManager() {};

    // Convert latitude and longitude to mesh code.
    static std::string LatLonToMeshCode(double lat, double lon);

    // Convert mesh code to latitude and longitude.
    static std::pair<double, double> MeshCodeToLatLon(const std::string mesh_code);

    // Convert latitude and longitude to mesh code.
    static inline double GetLatDifference(double distance, double ref_lat);
    static inline double GetLonDifference(double distance, double ref_lat);

    // Calculate distance between two points.
    static inline double Distance(double lat1, double lon1, double lat2, double lon2);
    static inline double MeshDistance(std::string mesh1, std::string mesh2);

    // Convert index to mesh code.
    static inline std::string IndexToMeshCode(int i, int j, int k, int nx, int ny, int nz, double dx, double dy, double dz, double lat0, double lon0);

    // Convert latitude and longitude to index.
    static inline std::tuple<int, int, int> LatLonDepthToIndex(double lat, double lon, double targetDepth, int nx, int ny, int nz, double dx, double dy, double dz, double lat0, double lon0, double depth0);

    // Convert index to latitude and longitude and depth.
    static inline std::tuple<double, double, double> IndexToLatLonDepth(int i, int j, int k, int nx, int ny, int nz, double dx, double dy, double dz, double lat0, double lon0, double depth0);

};

inline std::string CoordinateManager::LatLonToMeshCode(double lat, double lon)
{
    // int(mt.floor(lat*240       / 160)) * 10000000 + int(mt.floor((lon-100)*160       / 160)) * 100000 +    # 1次メッシュ
    //             int(mt.floor(lat*240 % 160 / 20))  * 10000    + int(mt.floor((lon-100)*160 % 160 / 20))  * 1000   +    # 2次メッシュ
    //             int(mt.floor(lat*240 % 20  / 2))   * 100      + int(mt.floor((lon-100)*160 % 20  / 2))   * 10     +    # 3次メッシュ
    //             int(mt.floor(lat*240)) % 2         * 2        + int(mt.floor((lon-100)*160)) % 2                  + 1) # 4次メッシュ
    
    int first_mesh = static_cast<int>(std::floor(lat * 240 / 160)) * 10000000 + static_cast<int>(std::floor((lon - 100) * 160 / 160)) * 100000;
    int second_mesh = static_cast<int>(std::floor(fmod(lat * 240, 160) / 20)) * 10000 + static_cast<int>(std::floor(fmod((lon - 100) * 160, 160) / 20)) * 1000;
    int third_mesh = static_cast<int>(std::floor(fmod(lat * 240, 20) / 2)) * 100 + static_cast<int>(std::floor(fmod((lon - 100) * 160, 20) / 2)) * 10;
    int fourth_mesh = static_cast<int>(std::floor(lat * 240)) % 2 * 2 + static_cast<int>(std::floor((lon - 100) * 160)) % 2 + 1;

    int mesh_code = first_mesh + second_mesh + third_mesh + fourth_mesh;

    return std::to_string(mesh_code);
}

inline std::pair<double, double> CoordinateManager::MeshCodeToLatLon(const std::string mesh_code)
{
    // メッシュコードから緯度経度を計算

    if (mesh_code.length() != 8)
    {
        std::cerr << "Invalid mesh code." << std::endl;
        return std::make_pair(-1, -1);
    }

    int lat_code1 = std::stoi(mesh_code.substr(0, 2));
    int lat_code2 = std::stoi(mesh_code.substr(4, 1));
    int lat_code3 = std::stoi(mesh_code.substr(6, 1));

    int lon_code1 = std::stoi(mesh_code.substr(2, 2));
    int lon_code2 = std::stoi(mesh_code.substr(5, 1));
    int lon_code3 = std::stoi(mesh_code.substr(7, 1));

    double lat = (lat_code1 * 80.0 + lat_code2 * 10.0 + lat_code3) * 30.0 / 3600.0 + 30.0 / 3600.0 / 2.0;
    double lon = (lon_code1 * 80.0 + lon_code2 * 10.0 + lon_code3) * 45.0 / 3600.0 + 100.0 + 45.0 / 3600.0 / 2.0;

    return std::make_pair(lat, lon);
}

inline double CoordinateManager::GetLatDifference(double distance, double ref_lat)
{
    // self.lat_per_km = abs(1/ (Ry / 360.0))#緯度/km

    // 緯度の差を計算
    double lat_per_km = 1.0 / (Ry / 360.0);
    double lat_diff = lat_per_km * distance / 1000.0;
    // abs
    if (lat_diff < 0)
    {
        lat_diff = -lat_diff;
    }
    return lat_diff;
}

inline double CoordinateManager::GetLonDifference(double distance, double ref_lat)
{

    // self.lon_per_km = abs(1/ (Rx * np.cos(np.radians(self.epi_lat))/360.0)) #経度/km      #ここにlonを入れてた注意

    // 参照緯度の余弦を使用して経度の差を計算
    double lon_per_km = 1.0 / (Rx * std::cos(ref_lat * M_PI / 180.0) / 360.0);
    double lon_diff = lon_per_km * distance / 1000.0;
    // abs
    if (lon_diff < 0)
    {
        lon_diff = -lon_diff;
    }
    return lon_diff;
}

inline double CoordinateManager::Distance(double lat1, double lon1, double lat2, double lon2)
{
    double ra = 6378140;       // equatorial radius (m)
    double rb = 6356755;       // polar radius (m)
    double F = (ra - rb) / ra; // flattening of the earth
    double rad_lat_a = lat1 * M_PI / 180.0;
    double rad_lon_a = lon1 * M_PI / 180.0;
    double rad_lat_b = lat2 * M_PI / 180.0;
    double rad_lon_b = lon2 * M_PI / 180.0;
    double pa = std::atan(rb / ra * std::tan(rad_lat_a));
    double pb = std::atan(rb / ra * std::tan(rad_lat_b));
    double xx = std::acos(std::sin(pa) * std::sin(pb) + std::cos(pa) * std::cos(pb) * std::cos(rad_lon_a - rad_lon_b));
    double c1 = (std::sin(xx) - xx) * std::pow(std::sin(pa) + std::sin(pb), 2) / std::pow(std::cos(xx / 2), 2);
    double c2 = (std::sin(xx) + xx) * std::pow(std::sin(pa) - std::sin(pb), 2) / std::pow(std::sin(xx / 2), 2);
    double dr = F / 8 * (c1 - c2);
    double rho = ra * (xx + dr);
    return rho;
}

inline double CoordinateManager::MeshDistance(std::string mesh1, std::string mesh2)
{
    // メッシュコードから緯度経度を計算
    std::pair<double, double> latlon1 = MeshCodeToLatLon(mesh1);
    std::pair<double, double> latlon2 = MeshCodeToLatLon(mesh2);

    // 緯度経度から距離を計算
    double distance = CoordinateManager::Distance(latlon1.first, latlon1.second, latlon2.first, latlon2.second);

    return distance;
}

inline std::string CoordinateManager::IndexToMeshCode(int i, int j, int k, int nx, int ny, int nz, double dx, double dy, double dz, double lat0, double lon0)
{
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

    double targetLat = lat0 - (maxLat - minLat) / 2.0 + static_cast<double>(i) / nx * deltaLat;
    double targetLon = lon0 - (maxLon - minLon) / 2.0 + static_cast<double>(j) / ny * deltaLon;

    // i,jに対応するメッシュコードを計算します。
    std::string targetMeshCode = CoordinateManager::LatLonToMeshCode(targetLat, targetLon);

    return targetMeshCode;
}

inline std::tuple<int, int, int> CoordinateManager::LatLonDepthToIndex(double targetLat, double targetLon, double targetDepth, int nx, int ny, int nz, double dx, double dy, double dz, double lat0, double lon0, double depth0)
{
    // 最も近いグリッドを探索します。
    // Argument Description
    // lat, lon, targetDepth: 緯度経度深度
    // nx, ny, nz: メッシュ数
    // dx, dy, dz: メッシュサイズ
    // lat0, lon0, depth0: メッシュの原点の緯度経度深度
    
    double minDistance = std::numeric_limits<double>::max();// double型の最大値
    int ii = -99, jj = -99, kk = -99;

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

    for (int i = 0; i < nx; ++i)
    {
        double gridLat = lat0 - (maxLat - minLat) / 2.0 + static_cast<double>(i) / nx * deltaLat;

        for (int j = 0; j < ny; ++j)
        {    
            double gridLon = lon0 - (maxLon - minLon) / 2.0 + static_cast<double>(j) / ny * deltaLon;

            double distanceHorizon = CoordinateManager::Distance(targetLat, targetLon, gridLat, gridLon);

            for (int k = 0; k < nz; ++k)
            {
                double gridDepth = depth0 + k * dz;    
                double distance = std::sqrt(distanceHorizon * distanceHorizon + (gridDepth - targetDepth) * (gridDepth - targetDepth));

                if (distance < minDistance)
                {
                    minDistance = distance;
                    ii = i;
                    jj = j;
                    kk = k;
                }
            }
        }
    }
    return std::make_tuple(ii, jj, kk);
}

inline std::tuple<double, double, double> CoordinateManager::IndexToLatLonDepth(int i, int j, int k, int nx, int ny, int nz, double dx, double dy, double dz, double lat0, double lon0, double depth0)
{
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

    double targetLat = lat0 - (maxLat - minLat) / 2.0 + static_cast<double>(i) / nx * deltaLat;
    double targetLon = lon0 - (maxLon - minLon) / 2.0 + static_cast<double>(j) / ny * deltaLon;
    double targetDepth = depth0 + k * dz;

    return std::make_tuple(targetLat, targetLon, targetDepth);
}
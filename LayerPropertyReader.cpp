#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <omp.h>
#include "CoordinateManager.h"
#include "LayerPropertyReader.h"

LayerPropertyReader::LayerPropertyReader(const std::string &filename, double minLat, double maxLat, double minLon, double maxLon)
    : file(filename), minLat(minLat), maxLat(maxLat), minLon(minLon), maxLon(maxLon)
{
    loadData();
}

void LayerPropertyReader::loadData()
{
    std::ifstream inputStream(file);
    if (!inputStream.is_open())
    {
        throw std::runtime_error("Unable to open file: " + file);
    }

    std::string line;
    while (std::getline(inputStream, line))
    {
        if (line[0] == '#')
            continue; // Ignore comment lines

        LayerProperty row = parseLine(line);
        std::pair<double, double> latlon = CoordinateManager::MeshCodeToLatLon(row.CODE);

        if (latlon.first >= minLat && latlon.first <= maxLat && latlon.second >= minLon && latlon.second <= maxLon)
        {
            data.push_back(row);
        }
    }
}

LayerProperty LayerPropertyReader::parseLine(const std::string &line)
{
    std::stringstream lineStream(line);
    std::string cell;
    LayerProperty row;

    std::getline(lineStream, cell, ',');
    row.CODE = cell.substr(0, 8);

    std::getline(lineStream, cell, ',');
    row.S0 = std::stod(cell);

    for (int i = 0; i < 30; ++i)
    {
        std::getline(lineStream, cell, ',');
        row.E[i] = std::stod(cell);
    }

    return row;
}

LayerProperty LayerPropertyReader::findLayerPropertyByCode(const std::string &meshcode)
{
    for (const auto &row : data)
    {
        if (row.CODE == meshcode)
        {
            return row;
        }
    }
    return LayerProperty();
}

LayerProperty LayerPropertyReader::findNearestLayerProperty(const std::string &meshcode)
{
    double minDistance = std::numeric_limits<double>::max();
    std::pair<double, double> targetLatLon = CoordinateManager::MeshCodeToLatLon(meshcode);
    LayerProperty nearest;

#pragma omp parallel
    {
        LayerProperty local_nearest;
        double local_minDistance = std::numeric_limits<double>::max();

#pragma omp for
        for (size_t i = 0; i < data.size(); ++i)
        {
            const auto &row = data[i];
            std::pair<double, double> latlon = CoordinateManager::MeshCodeToLatLon(row.CODE);
            double distance = CoordinateManager::Distance(latlon.first, latlon.second, targetLatLon.first, targetLatLon.second);
            if (local_minDistance > distance)
            {
                local_minDistance = distance;
                local_nearest = row;
            }
        }

#pragma omp critical
        {
            if (minDistance > local_minDistance)
            {
                minDistance = local_minDistance;
                nearest = local_nearest;
            }
        }
    }
    return nearest;
}

LayerProperty LayerPropertyReader::operator[](std::string meshcode)
{
    LayerProperty result = findLayerPropertyByCode(meshcode);

    if (result.CODE.empty())
    {
        result = findNearestLayerProperty(meshcode);
    }

    return result;
}

int LayerPropertyReader::getLayerByElevation(LayerProperty targetMesh, double targetElevation) const
{
    // LayerPropertyは上向が正、最も深い層（elevationが最も小さい値の層）を求める

    int STN = 33; // 最深層

    // EO が最も浅い層
    // E30 が最も深い層
    const int numberOfLayers = 31;
    double elevations[numberOfLayers];
    for (int i = 0; i < numberOfLayers; ++i)
    {
        elevations[i] = targetMesh.E[i];
    }

    // std::cout << "targetElevation = " << targetElevation << std::endl;

    // targetElevation
    for (int i = 0; i < numberOfLayers; ++i)
    {
        // std::cout << "elevations[" << i << "] = " << elevations[i] << std::endl;

        if (targetElevation > elevations[i])
        {
            STN = i;
            break;
        }
    }

    // if STN == 0, targetElevation is higher than the highest elevation
    if (STN == 0)
    {
        STN = 1;
        // std::cout << "targetElevation is higher than the highest elevation" << std::endl;
        // std::cout << "targetElevation = " << targetElevation << std::endl;
        // std::cout << "targetMesh.CODE = " << targetMesh.CODE << std::endl;
    }

    return STN;
}

// Receiver.cpp
#include "Receiver.h"
#include <fstream>
#include <iostream>
#include <string>
Receiver::Receiver(int id, double lat, double lon, double depth, std::string name, int x, int y, int z) : 
                    id(id), lat(lat), lon(lon), depth(depth), name(name),
                    x(x), y(y), z(z) {}

void Receiver::setReceiverParameters(int x, int y, int z)
{
    this->x = x;
    this->y = y;
    this->z = z;
}

void Receiver::record(double ux, double uy, double uz)
{
    // 受信点での速度を記録
    uxData.push_back(ux);
    uyData.push_back(uy);
    uzData.push_back(uz);
}

void Receiver::getLocation(int &x, int &y, int &z) const
{
    x = this->x;
    y = this->y;
    z = this->z;
}

void Receiver::getVelocityData(std::vector<double> &uxData, std::vector<double> &uyData, std::vector<double> &uzData) const
{
    // 受信点での速度を取得
    uxData = this->uxData;
    uyData = this->uyData;
    uzData = this->uzData;
}
void Receiver::getReceiverInfo(int &id, double &lat, double &lon, double &depth, std::string &name) const
{
    id = this->id;
    lat = this->lat;
    lon = this->lon;
    depth = this->depth;
    name = this->name;
}
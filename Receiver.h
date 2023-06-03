// Receiver.h
#include <vector>
#include <string>
#ifndef RECEIVER_H
#define RECEIVER_H

class Receiver {
public:
    Receiver(int id,double lat,double lon,double depth,std::string name, int x, int y, int z);

    void setReceiverParameters(int x, int y, int z);
    void record(double ux, double uy, double uz);
    void getLocation(int &x, int &y, int &z) const;
    void getVelocityData(std::vector<double> &uxData, std::vector<double> &uyData, std::vector<double> &uzData) const;
    void getReceiverInfo(int &id, double &lat, double &lon, double &depth, std::string &name) const;

private:
    int x, y, z;
    int id;
    double lat;
    double lon;
    double depth;// m
    std::string name;

    std::vector<double> uxData;
    std::vector<double> uyData;
    std::vector<double> uzData;
};


#endif // RECEIVER_H

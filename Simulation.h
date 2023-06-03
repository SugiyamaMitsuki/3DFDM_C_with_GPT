#ifndef SIMULATION_H
#define SIMULATION_H

#include "Grid.h"
#include "Source.h"
#include "Receiver.h"
#include <vector>

class Simulation
{
public:
    explicit Simulation(const std::string &configFile = "./input/config.csv");

    void initializeModel();
    void setPhysicalProperty(bool isTest, std::string inputPhysicalProperty, std::string inputLayerProperty);
    void savePhysicalProperty(std::string outputPhysicalProperty);
    void addSource(const Source &source);
    void setSource();
    void addReceiver(const Receiver &receiver);
    void setReceiver();
    void run();
    void writeData(const std::string &filename);

private:
    Grid grid;
    
    std::vector<Source> sources;
    std::vector<Receiver> receivers;

    double dt;
    int num_steps;
    double totalTime;
    int num_threads;
    int snapshotInterval;
    std::string snapshotFilename;

    int nx, ny, nz;
    double dx, dy, dz;
    double lat0, lon0, depth0;

    double f0;

    double spatialWeightNonReflective;
    double temporalWeightNonReflective;
    int absorbingGridNumber;
    double absorbingDampingFactor;

    void recordSnapshot(int time_step);
    void applySourceAtTimeStep(double dt, int time_step);
    void recordReceiversAtTimeStep(int time_step);
    void readConfig(const std::string &configFile);
    
};

#endif // SIMULATION_H

#include <iostream>

#include "Simulation.h"

int main()
{
  // This program is a 3D Finite Difference Method (FDM) simulation of earthquake wave simulation.
  // written by Mitsuki Sugiyama, 2023/04 with c++
  // This program is based on the following paper:
  // X = NS (North-South) direction (positive northward)
  // Y = EW (East-West) direction (positive eastward)
  // Z = UD (Up-Down) direction (positive downward)

  std::cout << "Start" << std::endl;

  // シミュレーションオブジェクトの生成
  std::string configFile = "./input/config.csv";
  Simulation simulation(configFile);
  simulation.run();

  std::cout << "Simulation finished." << std::endl;

  return 0;
}

#include "topology.h"
#include <iostream>
#include <fstream>
#include <sstream>

void Topology::Print() {
  std::cout << "Topology size: " << atomName.size() << '\n';
}

void Topology::Read(const Config &config) {
  std::ifstream fp(config.TrajName());
  std::string line;
  std::getline(fp, line);
  std::getline(fp, line);

  std::string atomname;
  for (size_t i = 0; i != config.Size(); ++i) {
    std::getline(fp, line);
    std::istringstream iss(line);
    iss >> atomname;
    atomName.push_back(atomname);
  }
}

void Frame::Init(int natoms) {
  Step = -1;
  Coords.resize(natoms);
}

void Frame::Read(std::ifstream &fp) {
  std::string line;
  std::getline(fp, line);
  std::getline(fp, line);

  std::string dummy;
  for (auto &atom : Coords) {
    std::getline(fp, line);
    std::istringstream iss(line);
    iss >> dummy >> atom[0] >> atom[1] >> atom[2];
  }
  ++Step;
}

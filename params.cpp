#include "params.h"
#include "topology.h"

#include <fstream>
#include <iostream>
#include <ostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_map>

Config::Config(std::string filename) {
  std::ifstream ifile(filename);

  if (!ifile.is_open()) {
    throw std::runtime_error("Config file not found!");
  }

  std::string line;
  std::string key;
  std::unordered_map<std::string, std::string> kv;

  while (std::getline(ifile, line)) {
    if (line.empty() || line[0] == '#') continue;
    std::istringstream iss(line);
    iss >> key;
    std::string rest;
    std::getline(iss, rest);
    if (!rest.empty() && rest[0] == ' ') rest.erase(0,1);
    kv[key] = rest;
  }

  ifile.close();

  fname = kv.at("filename");
  natoms = std::stoi(kv.at("natoms"));
  stride = std::stoi(kv.at("stride"));
  dt = std::stof(kv.at("dt"));
  group1 = kv.at("group1");
  group2 = kv.at("group2");

  std::istringstream iss(kv.at("atomtypes"));
  std::string type;
  while(iss >> type) atomtypes.push_back(type);

  iss.clear();
  iss.str(kv.at("axis-a"));
  iss >> axis_a;
  iss.clear();
  iss.str(kv.at("axis-b"));
  iss >> axis_b;
  iss.clear();
  iss.str(kv.at("axis-c"));
  iss >> axis_c;
}

// void Density(const Frame &frame, const std::string &atomname, int nslabs) {
//   auto &atoms {frame.coords};
//   int slab;
//   std::vector<int> dens(nslabs);
//   for(size_t i = 0; i != atoms.size(); ++i) {
//     if (atomname == atoms[i].name) {
//       slab = static_cast<int>(atoms[i].rx[2]);
//       dens[slab]++;
//     }
//   }
//   for (size_t i = 0; i != dens.size(); ++i) {
//     std::cout << i << ' ' << dens[i] << '\n';
//   }
// }

void Select( std::vector<int> &list, const Frame &frame, float za, float zb) {

  list.clear();
  for (size_t i = 0; i != frame.Coords.size(); ++i) {

    if (frame.Coords[i][0] > za && frame.Coords[i][0] < zb) {
      list.push_back(i);
    }
  }

  // std::cout << "Range list: " << list.size() << "\n";
}

void Config::Print() {
  std::cout << "Config file\n\n";
  std::cout << "Filename: " << fname << '\n';

  std::cout << "Atom types: ";
  for (const auto &at : atomtypes) {
    std::cout << at << ' ';
  }
  std::cout << '\n';

  std::cout << "dt: " << dt << '\n';

  std::cout << "a-axis: " << axis_a << '\n';

  std::cout << "b-axis: " << axis_b << '\n';

  std::cout << "c-axis: " << axis_c << '\n';
}

std::vector<int> String2IntList(const std::string &str) {
  std::vector<int> list;

  std::stringstream ss(str);
  std::string token;
  while (ss >> token) {
    if (token.back() == ',') {
      token.pop_back();
    }

    if (token[0] == '-') {
      ss >> token;
      if (token.back() == ',') {
        token.pop_back();
      }
      int start = list.back() + 1;
      int end = std::stoi(token) - 1;
      for (int i = start; i <= end; ++i) {
        list.push_back(i);
      }
    } else {
      list.push_back(std::stoi(token) - 1);
    }
  }

  return list;
}



void Acceptors::Print() {
  std::cout << "Atoms in group: " << atoms.size() << '\n';
  for (auto var : atoms) {
    std::cout << var << '\n';
  }
}

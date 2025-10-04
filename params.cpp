#include "params.h"

#include <fstream>
#include <iostream>
#include <ostream>
#include <sstream>
#include <string>

Config::Config(std::string filename) {
  std::ifstream ifile;

  ifile.open(filename);
  if (!ifile.is_open()) {
    std::cerr << "File not open\n";
  }

  std::string line;

  std::getline(ifile, line);
  fname = line;

  std::getline(ifile, line);
  natoms = std::stoi(line);

  std::getline(ifile, line);
  ntypes = std::stoi(line);

  for (size_t i = 0; i != ntypes; ++i) {
    std::getline(ifile, line);
    atomtypes.push_back(line);
  }

  std::getline(ifile, line);
  group1 = line;

  std::getline(ifile, line);
  group2 = line;

  std::getline(ifile, line);
  stride = std::stoi(line);

  std::getline(ifile, line);
  dt = std::stof(line);

  std::getline(ifile, line);
  std::istringstream iss(line);
  iss >> box_a[0] >> box_a[1] >> box_a[2];

  std::getline(ifile, line);
  iss.clear();
  iss.str(line);
  iss >> box_b[0] >> box_b[1] >> box_b[2];

  std::getline(ifile, line);
  iss.clear();
  iss.str(line);
  iss >> box_c[0] >> box_c[1] >> box_c[2];

  std::getline(ifile, line);
  slabs = std::stoi(line);

  if (slabs == true) {
    std::getline(ifile, line);
    nslabs = std::stoi(line);
  }

  ifile.close();
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

void Config::Print() {
  std::cout << "Config file\n\n";
  std::cout << "Filename: " << fname << '\n';

  std::cout << "Atom types: ";
  for (const auto &at : atomtypes) {
    std::cout << at << ' ';
  }
  std::cout << '\n';

  std::cout << "dt: " << dt << '\n';

  std::cout << "a-axis: ";
  for (const auto &v : box_a) {
    std::cout << v << ' ';
  }
  std::cout << '\n';

  std::cout << "b-axis: ";
  for (const auto &v : box_b) {
    std::cout << v << ' ';
  }
  std::cout << '\n';

  std::cout << "c-axis: ";
  for (const auto &v : box_c) {
    std::cout << v << ' ';
  }
  std::cout << '\n';

  if (slabs == true) {
    std::cout << "nslabs: " << nslabs << '\n';
  }
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

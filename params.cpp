#include "params.h"
#include "topology.h"

#include <fstream>
#include <iostream>
#include <ostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_map>

Axis CharToAxis(char c) {
    switch (c) {
        case 'a': case 'A': case 'x': case 'X': return X;
        case 'b': case 'B': case 'y': case 'Y': return Y;
        case 'c': case 'C': case 'z': case 'Z': return Z;
        default:
            throw std::invalid_argument("Invalid axis character. Use a/b/c or x/y/z.");
    }
}

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

  auto parseKey = [&](const std::string &key, const std::string &def) -> std::string {
    auto it = kv.find(key);
    return (it != kv.end()) ? it->second : def;
  };

  fname = parseKey("filename", "XYZCAR");
  natoms = std::stoul(kv.at("natoms"));
  stride = std::stoul(parseKey("stride", "1"));
  dt = std::stof(kv.at("dt"));
  group1 = kv.at("group1");
  group2 = kv.at("group2");
  za = std::stof(kv.at("za"));
  zb = std::stof(kv.at("zb"));
  bframe = std::stoi(parseKey("bframe", "0"));
  eframe = std::stoi(parseKey("eframe", "-1"));
  std::string ax = kv.at("axis");
  axis = CharToAxis(ax[0]);

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

  std::cout << "range: ( " << za << " - " << zb << " )\n";

  std::cout << "Normal axis: " << axis << '\n';
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

#pragma once
#include "Math/Matrix.h"
#include <array>
#include <string>
#include <vector>

class Config {
private:
  std::string fname;
  int natoms;
  int ntypes;
  std::vector<std::string> atomtypes;
  std::string group1;
  std::string group2;
  int stride;
  float dt;
  float za;
  float zb;
  Math::Vec3 axis_a, axis_b, axis_c;

public:
  Config(std::string fname);
  void Print();
  std::string TrajName() const { return fname; }
  int Size() const { return natoms; }
  std::string GetGroup1() const { return group1; }
  std::string GetGroup2() const { return group2; }
  float getZa() const { return za; }
  float getZb() const { return zb; }
  Math::Matrix33 Box() const { return Math::Matrix33 {axis_a,axis_b,axis_c}; }
};

std::vector<int> String2IntList(const std::string &str);

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
  Math::Vec3 box_a, box_b, box_c;
  bool slabs;
  int nslabs;

public:
  Config(std::string fname);
  void Print();
  std::string TrajName() const { return fname; }
  int Size() const { return natoms; }
  std::string GetGroup1() const { return group1; }
  std::string GetGroup2() const { return group2; }
  Math::Matrix33 Box() const { return Math::Matrix33 {box_a,box_b,box_c}; }
};

class Acceptors {
  // static std::vector<std::string> acc;
  std::vector<int> atoms;

  void Print();
};

std::vector<int> ParseList(const std::vector<int> &list);

std::vector<int> String2IntList(const std::string &str);

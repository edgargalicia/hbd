#pragma once
#include "Math/Matrix.h"
#include <array>
#include <string>
#include <vector>

enum Axis { X = 0, Y = 1, Z = 2 };

class Config {
private:
  std::string fname;
  size_t natoms;
  std::vector<std::string> atomtypes;
  std::string group1;
  std::string group2;
  size_t stride;
  float dt;
  float za;
  float zb;
  Math::Vec3 axis_a, axis_b, axis_c;
  int bframe;
  int eframe;
  Axis axis;

public:
  Config(std::string fname);
  void Print();
  std::string TrajName() const { return fname; }
  size_t Size() const { return natoms; }
  std::string GetGroup1() const { return group1; }
  std::string GetGroup2() const { return group2; }
  float getZa() const { return za; }
  float getZb() const { return zb; }
  Math::Matrix33 Box() const { return Math::Matrix33 {axis_a,axis_b,axis_c}; }
  int Begin() const { return bframe; }
  int End() const { return eframe; }
  void SeteFrames(int e) { eframe = e; }
  Axis GetAxis() const { return axis; }
};

std::vector<int> String2IntList(const std::string &str);

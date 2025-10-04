#include "params.h"

#include <vector>
#include <string>

struct Topology {
private:
  std::vector<std::string> atomName;

public:
  void Print();
  void Read(const Config &config);
};

struct Frame {
  int Step;
  std::vector<std::array<float, 3>> Coords;

  void Init(int natoms);
  void Read(std::ifstream &fp);
};

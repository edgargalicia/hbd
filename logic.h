#include "params.h"
#include "hbmap.h"
#include "topology.h"
#include <string>
#include <fstream>

class Logic {
private:
  std::string filename;
  Config conf;
  Topology topo;
  Frame frame;
  Box box;
  HBMap hbmap;
  std::ifstream infile;

public:
  Logic() : filename{"input.dat"}, conf{filename} {}
  void Init();
  void Run();
  void Finish();
};


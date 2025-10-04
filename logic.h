#include "params.h"
#include "topology.h"
#include <string>
#include <fstream>

class Logic {
private:
  std::string filename;
  Config conf;
  Topology topo;
  Frame frame;
  std::ifstream infile;

public:
  Logic();
  void Init();
  void Run();
  void Finish();
};


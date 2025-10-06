#include "params.h"

#include <Math/Matrix.h>
#include <vector>
#include <string>

struct Box {
  Math::Matrix33 mbox;
  Math::Vec3 fbox;
  Math::Vec3 hbox;
  Math::Vec3 mhbox;

  void Pbc( Math::Vec3 &dx ) const;
};

Box InitBox( const Math::Matrix33 &mbox );

class Topology {
private:
  std::vector<std::string> atomNames;
  std::vector<std::pair<int, int>> bonds;

public:
  void Read(const Config &config, const Box &box);
  void Print();
  std::vector<std::string> AtomName() const { return atomNames; }
  std::vector<std::pair<int, int>> Bonds() const { return bonds; }
};

struct Frame {
  int Step;
  std::vector<Math::Vec3> Coords;

  void Init(int natoms);
  bool Read(std::ifstream &fp);
};

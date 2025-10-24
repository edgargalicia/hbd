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
  // void Read(const Config &config, const Box &box);
  void Print();
  void PushAtom(const std::string &at) { atomNames.push_back(at); }
  void PushBond(int i, int j) { bonds.push_back(std::make_pair(i, j)); }
  std::vector<std::string> AtomName() const { return atomNames; }
  std::vector<std::pair<int, int>> Bonds() const { return bonds; }
};

Topology ReadTopology(const Config &config, const Box &box);

class Frame {
private:
  int step;
  std::vector<Math::Vec3> coords;

public:
  void Init(int natoms);
  bool Read(std::ifstream &fp);
  size_t Size() const { return coords.size(); }
  int Step() const { return step; }
  const std::vector<Math::Vec3> &Coords() const { return coords; }
};

void Select( std::vector<int> &list, const Frame &frame, float za , float zb, size_t axis );

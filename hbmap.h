#include <params.h>
#include <string>
#include <unordered_map>
#include <vector>

class Topology;
class HBMap {
private:
  std::vector<int> acceptors;
  std::vector<int> mapAcceptors;
  std::unordered_map<int, std::vector<int>> donors;
  std::vector<int> donorOrder;
  std::vector<int> mapDonors;

public:
  HBMap() = default;
  HBMap(const std::string &str, const Topology &topo);
  ~HBMap() = default;
  int TotalProtons() const;
  void Print();
  void FindDonors(const Topology &topo);
  std::vector<int> getAcceptors() const { return acceptors; }
  std::vector<int> getDonors() const { return donorOrder; }
  std::vector<int> getAccMap() const { return mapAcceptors; }
  std::vector<int> getDonMap() const { return mapDonors; }
  const std::unordered_map<int, std::vector<int>>& getDonUno() const { return donors; }
};

std::vector<int> Mapping(const std::vector<int> &list, int natoms);
class Topology;
std::vector<int> Match(const std::vector<int> &atomList, const Topology &topology, const std::string &str);

class Vec3;
class Box;
int isHBonded( int d, int h, int a, std::vector<Math::Vec3> &x, const Box &box, float &d_ha, float &ang );

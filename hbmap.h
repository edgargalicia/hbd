#include <params.h>
#include <string>
#include <unordered_map>
#include <vector>

class Topology;
class HBMap {
private:
  std::vector<int> Acceptors;
  std::vector<int> MapAcceptors;
  std::unordered_map<int, std::vector<int>> Donors;
  std::vector<int> MapDonors;

public:
  HBMap(const std::string &str, const Topology &topo);
  ~HBMap() = default;
  int TotalProtons() const;
  void Print();
  void FindDonors(const Topology &topo);
  std::unordered_map<int, std::vector<int>> DonorToProtons() const { return Donors; }
};

std::vector<int> Mapping(const std::vector<int> &list, int natoms);
class Topology;
std::vector<int> Match(const std::vector<int> &atomList, const Topology &topology, const std::string &str);

class Vec3;
class Box;
bool isHBonded( int d, int h, int a, Vec3 *x[], const Box &box, float &d_ha, float &ang );

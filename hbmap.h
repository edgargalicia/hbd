#include <fstream>
#include <functional>
#include <params.h>
#include <string>
#include <unordered_map>
#include <unordered_set>
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
  HBMap(const std::vector<int> &list, const Topology &topo);
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
bool isHBonded( int don, int prot, int acc, std::vector<Math::Vec3> &x, const Box &box, float &d_ha, float &ang );
void PrintTimeSeries( std::ofstream &fp, int acc, int don, int prot, float d_ha, float ang );
void PrintStats( std::ofstream &fp, int hbond, int intra1, int intra2, int inter, const std::unordered_set<int> &activeAcceptors, const std::unordered_set<int> &activeDonors, const std::unordered_set<int> &group1Set, const std::unordered_set<int> &group2Set );

struct BondKey {
  int acceptor;
  int donor;

  bool operator==(const BondKey &other) const {
    return acceptor == other.acceptor && donor == other.donor;
  }
};

struct BondKeyHash {
  std::size_t operator()(const BondKey &key) const noexcept {
    return std::hash<int>()(key.acceptor) ^ (std::hash<int>()(key.donor) << 1);
  }
};

enum class BondType { Intra1, Intra2, Inter};

struct BondInfo {
  BondKey key;
  BondType type;
};

void addBondPresence(std::unordered_map<BondKey, std::vector<int>, BondKeyHash> &bondPresence,
                     std::unordered_map<BondKey, BondType, BondKeyHash> &bondTypeMap,
                     int acc, int don, int frame,
                     const std::unordered_set<int> &group1Set,
                     const std::unordered_set<int> &group2Set);

std::unordered_map<BondKey, double, BondKeyHash> computeBondPresencePercentages(
  const std::unordered_map<BondKey, std::vector<int>, BondKeyHash> &bonds, int totalFrames);

struct GroupActivity {
  int g1Don = 0, g1Acc = 0, g2Don = 0, g2Acc = 0;
};

GroupActivity CountGroupActivity(const std::unordered_set<int> &activeAcceptors, const std::unordered_set<int> &activeDonors, const std::unordered_set<int> &group1Set, const std::unordered_set<int> &group2Set);

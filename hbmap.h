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

  void findDonors(const Topology &topo);

public:
  HBMap() = default;
  HBMap(const std::vector<int> &list, const Topology &topo);
  ~HBMap() = default;
  int TotalProtons() const;
  void Print();
  std::vector<int> getAcceptors() const { return acceptors; }
  std::vector<int> getDonors() const { return donorOrder; }
  std::vector<int> getAccMap() const { return mapAcceptors; }
  std::vector<int> getDonMap() const { return mapDonors; }
  const std::unordered_map<int, std::vector<int>>& getDonUno() const { return donors; }
};

struct BondPair {
  int acceptor;
  int donor;

  bool operator==(const BondPair &other) const {
    return acceptor == other.acceptor && donor == other.donor;
  }
};

struct BondKeyHash {
  std::size_t operator()(const BondPair &key) const noexcept {
    return std::hash<int>()(key.acceptor) ^ (std::hash<int>()(key.donor) << 1);
  }
};

enum class BondType { Intra1, Intra2, Inter};

struct BondInfo {
  BondPair key;
  BondType type;
};

std::unordered_map<BondPair, double, BondKeyHash> computeBondPresencePercentages(
  const std::unordered_map<BondPair, std::vector<int>, BondKeyHash> &bonds, int totalFrames);

struct GroupActivity {
  int g1Don = 0, g1Acc = 0, g2Don = 0, g2Acc = 0;
};

GroupActivity CountGroupActivity(const std::unordered_set<int> &activeAcceptors, const std::unordered_set<int> &activeDonors, const std::unordered_set<int> &group1Set, const std::unordered_set<int> &group2Set);

struct BondAtoms {
  BondPair bondpair;
  int proton;
};

struct HBGeometry {
    float d_da;  // Hydrogen–Acceptor distance
    float d_ha;  // Hydrogen–Acceptor distance
    float ang;   // Donor–Hydrogen–Acceptor angle
};

struct HBBinner {
  float d_min = 0.0f, d_max = 4.0f, d_bin = 0.1f;
  float a_min = 0.0f, a_max = 90.0f, a_bin = 1.0f;
  int n_d = int((d_max - d_min) / d_bin);
  int n_a = int((a_max - a_min) / a_bin);
  std::vector<int> count_dha;
  std::vector<int> count_dda;
  std::vector<int> count_a;

  HBBinner() : count_dha(n_d,0), count_dda(n_d,0), count_a(n_a,0) {}

  void add(const HBGeometry &g) {
    int i = int((g.d_da - d_min) / d_bin);
    int j = int((g.d_ha - d_min) / d_bin);
    if (i >= 0 && i < n_d) {
      count_dda[i]++;
      count_dha[j]++;
    }
    i = int((g.ang  - a_min) / a_bin);
    if (i >= 0 && i < n_a)
      count_a[i]++;
  }

  void print(std::ostream &os) const {
    os << "# Distribution of distances D-A\tH-A\n";
    for (int i = 0; i < n_d; ++i) {
      float d_center = d_min + (i + 0.5f) * d_bin;
      os << d_center << "\t" << count_dda[i] << "\t" << count_dha[i] << '\n';
    }
    os << "\n\n";
    os << "# Distribution of angles\n";
    for (int i = 0; i < n_a; ++i) {
      float a_center = a_min + (i + 0.5f) * a_bin;
      os << a_center << "\t" << count_a[i] << '\n';
    }
  }
};

struct Box;
struct HBondContext {
    const std::vector<Math::Vec3> &coords;
    const Box &box;
};

struct HBondStats {
  int hbond;
  int intra1;
  int intra2;
  int inter;
  std::unordered_set<int> activeAcceptors;
  std::unordered_set<int> activeDonors;

  void reset();
};

bool isHBonded( const BondAtoms &at, const HBondContext &ctx, HBGeometry &hbg );
void PrintTimeSeries( std::ofstream &fp, int acc, int don, int prot, float d_ha, float ang );
void PrintStats( std::ofstream &fp, const HBondStats &hbs, const std::unordered_set<int> &group1Set, const std::unordered_set<int> &group2Set);

void addBondPresence(std::unordered_map<BondPair, std::vector<int>, BondKeyHash> &bondPresence,
                     std::unordered_map<BondPair, BondType, BondKeyHash> &bondTypeMap,
                     const BondAtoms &at, int frame,
                     const std::unordered_set<int> &group1Set,
                     const std::unordered_set<int> &group2Set);

